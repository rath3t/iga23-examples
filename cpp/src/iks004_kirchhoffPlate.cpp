// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <matplot/matplot.h>
#include <numbers>

#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/scalarFE.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/duneUtilities.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/loadControlObserver.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

using namespace Ikarus;
template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
class KirchhoffPlate : public ScalarFieldFE<typename Basis_::FlatBasis> {
 public:
  using Basis             = Basis_;
  using FlatBasis         = typename Basis::FlatBasis;
  using BaseDisp          = ScalarFieldFE<FlatBasis>;  // Handles globalIndices function
  using LocalView         = typename FlatBasis::LocalView;
  using Element           = typename LocalView::Element;
  using Geometry          = typename Element::Geometry;
  using FERequirementType = FERequirements_;
  using Traits            = TraitsFromLocalView<LocalView, useEigenRef>;

  KirchhoffPlate(const Basis &basis, const typename LocalView::Element &element, double p_Emodul, double p_nu,
                 double p_thickness)
      : BaseDisp(basis.flat(), element), Emodul{p_Emodul}, nu{p_nu}, thickness{p_thickness} {
    this->localView().bind(element);
    geometry_.emplace(this->localView().element().geometry());
  }

  static Eigen::Matrix<double, 3, 3> constitutiveMatrix(double Emod, double p_nu, double p_thickness) {
    const double factor = Emod * Dune::power(p_thickness, 3) / (12.0 * (1.0 - p_nu * p_nu));
    Eigen::Matrix<double, 3, 3> D;
    D.setZero();
    D(0, 0) = 1;
    D(0, 1) = D(1, 0) = p_nu;
    D(1, 1)           = 1;
    D(2, 2)           = (1 - p_nu) / 2.0;
    D *= factor;
    return D;
  }

  inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }

 protected:
  template <typename ScalarType>
  auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
                                                         = std::nullopt) const -> ScalarType {
    const auto &wGlobal = par.getGlobalSolution(Ikarus::FESolutions::displacement);
    const auto &lambda  = par.getParameter(Ikarus::FEParameter::loadfactor);
    const auto D        = constitutiveMatrix(Emodul, nu, thickness);
    ScalarType energy   = 0.0;
    auto &ele           = this->localView().element();
    auto &fe            = this->localView().tree().finiteElement();
    Eigen::VectorX<ScalarType> wNodal;
    wNodal.setZero(fe.size());
    Dune::CachedLocalBasis localBasis(fe.localBasis());
    const auto &rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * localBasis.order());

    localBasis.bind(rule, Dune::bindDerivatives(0, 2));
    if (dx) {
      for (auto i = 0U; i < fe.size(); ++i)
        wNodal(i) = dx.value()[i] + wGlobal[this->localView().index(this->localView().tree().localIndex(i))[0]];
    } else {
      for (auto i = 0U; i < fe.size(); ++i)
        wNodal(i) = wGlobal[this->localView().index(this->localView().tree().localIndex(i))[0]];
    }

    /// Calculate Kirchhoff plate energy
    for (auto &&[gpIndex, gp] : localBasis.viewOverIntegrationPoints()) {
      auto &N          = localBasis.evaluateFunction(gpIndex);
      auto &ddN        = localBasis.evaluateSecondDerivatives(gpIndex);
      auto &ddN_xixi   = ddN.col(0);
      auto &ddN_etaeta = ddN.col(1);
      auto &ddN_xieta  = ddN.col(2);

      const auto Jinv = Dune::toEigen(geometry_->jacobianInverseTransposed(gp.position())).transpose().eval();

      Eigen::VectorXd ddN_xx(fe.size());
      Eigen::VectorXd ddN_yy(fe.size());
      Eigen::VectorXd ddN_xy(fe.size());
      using Dune::power;
      // The following derivative transformation assumes a non-distorted grid, otherwise there would be non-linear terms
      for (auto i = 0U; i < fe.size(); ++i) {
        ddN_xx[i] = ddN_xixi[i] * power(Jinv(0, 0), 2);
        ddN_yy[i] = ddN_etaeta[i] * power(Jinv(1, 1), 2);
        ddN_xy[i] = ddN_xieta[i] * Jinv(0, 0) * Jinv(1, 1);
      }
      Eigen::Vector<ScalarType, 3> kappa;
      kappa << ddN_xx.dot(wNodal), ddN_yy.dot(wNodal), 2 * ddN_xy.dot(wNodal);
      ScalarType w = N.dot(wNodal);

      energy += (0.5 * kappa.dot(D * kappa) - w * lambda) * geometry_->integrationElement(gp.position()) * gp.weight();
    }

    /// Clamp boundary using penalty method
    const double penaltyFactor = 1e8;
    if (ele.hasBoundaryIntersections())
      for (auto &intersection : intersections(this->localView().globalBasis().gridView(), ele))
        if (intersection.boundary()) {
          const auto &rule1 = Dune::QuadratureRules<double, 1>::rule(intersection.type(), 2 * localBasis.order());
          Eigen::MatrixX2d dN_xi_eta;
          for (auto &gp : rule1) {
            const auto &gpInElement = intersection.geometryInInside().global(gp.position());
            localBasis.evaluateJacobian(gpInElement, dN_xi_eta);
            Eigen::VectorXd dN_x(fe.size());
            Eigen::VectorXd dN_y(fe.size());
            const auto Jinv = Dune::toEigen(geometry_->jacobianInverseTransposed(gpInElement)).transpose().eval();
            for (auto i = 0U; i < fe.size(); ++i) {
              dN_x[i] = dN_xi_eta(i, 0) * Jinv(0, 0);
              dN_y[i] = dN_xi_eta(i, 1) * Jinv(1, 1);
            }
            const ScalarType w_x = dN_x.dot(wNodal);
            const ScalarType w_y = dN_y.dot(wNodal);

            energy += 0.0 * 0.5 * penaltyFactor * (w_x * w_x + w_y * w_y);
          }
        }

    return energy;
  }

 private:
  /// Dune::Geometry<...> is not copy assignable, see https://gitlab.dune-project.org/core/dune-grid/-/issues/140,
  /// Thus, we wrap it inside a std::optional
  std::optional<typename LocalView::Element::Geometry> geometry_;

  double Emodul;
  double nu;
  double thickness;
};

int main(int argc, char **argv) {
  Ikarus::init(argc, argv);

  /// Create 2D nurbs grid
  using namespace Ikarus;
  constexpr int griddim                                    = 2;
  constexpr int dimworld                                   = 2;
  const std::array<std::vector<double>, griddim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;

  const double Lx = 1;
  const double Ly = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {0, Ly}, .w = 1}}, {{.p = {Lx, 0}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, griddim> dimsize = {2, 2};

  std::vector<double> dofsVec;
  std::vector<double> l2Evcector;
  auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;

  Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  /// Increate polynomial degree in each direction
  patchData = Dune::IGA::degreeElevate(patchData, 0, 1);
  patchData = Dune::IGA::degreeElevate(patchData, 1, 1);
  Grid grid(patchData);

  for (int ref = 0; ref < 5; ++ref) {
    auto gridView = grid.leafGridView();
    //    draw(gridView);
    using namespace Dune::Functions::BasisFactory;
    /// Create nurbs basis with extracted preBase from grid
    auto basis = Ikarus::makeBasis(gridView, nurbs());
    /// Fix complete boundary (simply supported plate)
    auto basisP = std::make_shared<const decltype(basis)>(basis);
    Ikarus::DirichletValues dirichletValues(basisP->flat());
    dirichletValues.fixBoundaryDOFs(
        [&](auto &dirichletFlags, auto &&globalIndex) { dirichletFlags[globalIndex] = true; });

    /// Create finite elements
    auto localView         = basis.flat().localView();
    const double Emod      = 2.1e8;
    const double nu        = 0.3;
    const double thickness = 0.1;
    std::vector<AutoDiffFE<KirchhoffPlate<decltype(basis)>>> fes;
    for (auto &ele : elements(gridView))
      fes.emplace_back(basis, ele, Emod, nu, thickness);

    /// Create assembler
    auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

    /// Create non-linear operator with potential energy
    Eigen::VectorXd w;
    w.setZero(basis.flat().size());

    double totalLoad = 2000;

    auto req = FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);
    
    auto kFunction = [&](auto &&disp_, auto &&lambdaLocal) -> auto & {
      req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
          .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
      return sparseAssembler.getMatrix(req);
    };

    auto rFunction = [&](auto &&disp_, auto &&lambdaLocal) -> auto & {
      req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
          .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
      return sparseAssembler.getVector(req);
    };

    const auto &K = kFunction(w, totalLoad);
    const auto &R = rFunction(w, totalLoad);
    auto solver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);

    solver.compute(K);
    solver.solve(w,-R);

    // Output solution to vtk


    /// Create analytical solution function for the simply supported case
    const double D = Emod * Dune::power(thickness, 3) / (12 * (1 - Dune::power(nu, 2)));
    // https://en.wikipedia.org/wiki/Bending_of_plates#Simply-supported_plate_with_uniformly-distributed_load
    auto wAna = [&](auto x) {
      double w                = 0.0;
      const int seriesFactors = 40;
      const double pi         = std::numbers::pi;
      auto oddFactors
          = std::ranges::iota_view(1, seriesFactors) | std::views::filter([](auto i) { return i % 2 != 0; });
      for (auto m : oddFactors)
        for (auto n : oddFactors)
          w += sin(m * pi * x[0] / Lx) * sin(n * pi * x[1] / Ly)
               / (m * n * Dune::power(m * m / (Lx * Lx) + n * n / (Ly * Ly), 2));

      return 16 * totalLoad / (Dune::power(pi, 6) * D) * w;
    };
    //    std::cout << wxy(Lx / 2.0, Ly / 2.0) << std::endl;


    //    std::cout << wCenterClamped << std::endl;
    auto wGlobalAnalyticFunction = Dune::Functions::makeAnalyticGridViewFunction(wAna, gridView);
    auto wGlobalFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis.flat(), w);
    auto localw                  = localFunction(wGlobalFunction);
    auto localwAna               = localFunction(wGlobalAnalyticFunction);

    Dune::SubsamplingVTKWriter vtkWriter(gridView, Dune::refinementLevels(2));
    vtkWriter.addVertexData(wGlobalFunction, Dune::VTK::FieldInfo("wFE", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.addVertexData(wGlobalAnalyticFunction, Dune::VTK::FieldInfo("wANA", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("iks004_kirchhoffPlate_"+ std::to_string(ref));

        /// Displacement at center of clamped square plate
    // clamped sol http://faculty.ce.berkeley.edu/rlt/reports/clamp.pdf
    const double wCenterClamped = 1.265319087 / (D / (totalLoad * Dune::power(Lx, 4)) * 1000.0);

    /// Calculate L_2 error for simply supported case
    double l2_error = 0.0;
    for (auto &ele : elements(gridView)) {
      localView.bind(ele);
      localw.bind(ele);
      localwAna.bind(ele);
      const auto geo   = localView.element().geometry();
      const auto &rule = Dune::QuadratureRules<double, 2>::rule(
          ele.type(), 2 * localView.tree().finiteElement().localBasis().order());
      for (auto gp : rule) {
        const auto gpGlobalPos = geo.global(gp.position());

        const auto w_ex = localwAna(gp.position());
        const auto w_fe = localw(gp.position());
        l2_error += Dune::power(w_ex - w_fe, 2) * ele.geometry().integrationElement(gp.position()) * gp.weight();
      }
    }

    l2_error = std::sqrt(l2_error);
    std::cout << "l2_error: " << l2_error << " Dofs:: " << basis.flat().size() << std::endl;
    dofsVec.push_back(basis.flat().size());
    l2Evcector.push_back(l2_error);
    grid.globalRefine(1);
  }
  /// Draw L_2 error over dofs count
  using namespace matplot;
  auto f  = figure(true);
  auto ax = gca();
  ax->y_axis().label("L2_error");

  ax->x_axis().label("#Dofs");
  auto p = ax->loglog(dofsVec, l2Evcector);
  p->line_width(2);
  p->marker(line_spec::marker_style::asterisk);
  //  f->draw();
  save("ConvergenceplotKirchhoffPlate.png");
  using namespace std::chrono_literals;
  // std::this_thread::sleep_for(5s);
}
