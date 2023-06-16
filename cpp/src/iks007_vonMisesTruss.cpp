// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <matplot/matplot.h>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <autodiff/forward/dual/dual.hpp>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/powerBasisFE.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/finiteElements/physicsHelper.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/duneUtilities.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/genericControlObserver.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

using namespace Ikarus;
template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
class Truss : public PowerBasisFE<typename Basis_::FlatBasis> {
 public:
  using Basis             = Basis_;
  using FlatBasis         = typename Basis::FlatBasis;
  using BaseDisp          = Ikarus::PowerBasisFE<FlatBasis>;
  using LocalView         = typename FlatBasis::LocalView;
  using Element           = typename LocalView::Element;
  using Geometry          = typename Element::Geometry;
  using FERequirementType = FERequirements_;
  using Traits            = TraitsFromLocalView<LocalView, useEigenRef>;
  Truss(const Basis &basis, const typename LocalView::Element &element, double p_EA)
      : BaseDisp(basis.flat(), element), EA{p_EA} {
    this->localView().bind(element);
  }

  inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }

 protected:
  template <typename ScalarType>
  auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
                                                         = std::nullopt) const -> ScalarType {
    const auto &d      = par.getGlobalSolution(Ikarus::FESolutions::displacement);
    const auto &lambda = par.getParameter(FEParameter::loadfactor);

    auto &ele     = this->localView().element();
    const auto X1 = Dune::toEigen(ele.geometry().corner(0));
    const auto X2 = Dune::toEigen(ele.geometry().corner(1));

    Eigen::Matrix<ScalarType, Traits::worlddim, 2> u;
    u.setZero();
    if (dx) {
      for (int i = 0; i < 2; ++i)
        for (int k2 = 0; k2 < Traits::worlddim; ++k2)
          u.col(i)(k2) = dx.value()[Traits::worlddim * i + k2]
                         + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
    } else {
      for (int i = 0; i < 2; ++i)
        for (int k2 = 0; k2 < Traits::worlddim; ++k2)
          u.col(i)(k2) = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
    }

    const Eigen::Vector2<ScalarType> x1 = X1 + u.col(0);
    const Eigen::Vector2<ScalarType> x2 = X2 + u.col(1);

    const double LRefsquared  = (X1 - X2).squaredNorm();
    const ScalarType lsquared = (x1 - x2).squaredNorm();

    const ScalarType Egl = 0.5 * (lsquared - LRefsquared) / LRefsquared;

    return 0.5 * EA / sqrt(LRefsquared) * Egl * Egl;
  }

 private:
  double EA;
};

int main(int argc, char **argv) {
  Ikarus::init(argc, argv);
  /// Construct grid
  Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
  const double h = 1.0;
  const double L = 1.0;
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L, h});
  gridFactory.insertVertex({2 * L, 0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();
  // draw(gridView);

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  /// Create finite elements
  const double EA = 100;
  std::vector<AutoDiffFE<Truss<decltype(basis)>>> fes;
  for (auto &ele : elements(gridView))
    fes.emplace_back(basis, ele, EA);

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto &dirichletFlags, auto &&globalIndex) { dirichletFlags[globalIndex] = true; });

  /// Create assembler
  auto denseFlatAssembler = DenseFlatAssembler(fes, dirichletValues);

  /// Create non-linear operator
  double lambda = 0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto RFunction = [&](auto &&u, auto &&lambdaLocal) -> auto & {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, u)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    auto &R = denseFlatAssembler.getVector(req);
    R[3] -= -lambdaLocal;
    return R;
  };
  auto KFunction = [&](auto &&u, auto &&lambdaLocal) -> auto & {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, u)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return denseFlatAssembler.getMatrix(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(functions(RFunction, KFunction), parameter(d, lambda));

  /// Choose linear solver
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);

  /// Create Nonlinear solver for controlroutine, i.e. a Newton-Rahpson object
  auto nr = Ikarus::makeNewtonRaphson(nonLinOp, std::move(linSolver));
  nr->setup({.tol = 1e-8, .maxIter = 100});

  /// Create Observer to write information of the non-linear solver on the
  /// console
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  const int loadSteps = 10;
  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);
  /// Create Observer which executes when control routines messages
  /// SOLUTION_CHANGED
  auto lvkObserver = std::make_shared<Ikarus::GenericControlObserver>(ControlMessages::SOLUTION_CHANGED, [&](int step) {
    lambdaAndDisp(0, step) = lambda;
    lambdaAndDisp(1, step) = d[2];
    lambdaAndDisp(2, step) = d[3];
  });

  /// Create Observer which writes vtk files when control routines messages
  /// SOLUTION_CHANGED
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter->setFileNamePrefix("iks007_vonMisesTruss");
  nr->subscribeAll(nonLinearSolverObserver);

  /// Create loadcontrol
  auto lc = Ikarus::LoadControl(std::move(nr), loadSteps, {0, 30});
  lc.subscribeAll({vtkWriter, lvkObserver});

  /// Execute!
  lc.run();

  /// Postprocess
  using namespace matplot;
  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd dVec      = -lambdaAndDisp.row(2);
  auto f                    = figure(true);

  title("Load-Displacement Curve");
  xlabel("y-Displacement");
  ylabel("LoadFactor");

  auto analyticalLoadDisplacementCurve = [&](auto &w) {
    const double Ltruss = std::sqrt(h * h + L * L);
    return EA * Dune::power(h, 3) / Dune::power(Ltruss, 3)
           * (w / h - 1.5 * Dune::power(w / h, 2) + 0.5 * Dune::power(w / h, 3));
  };

  std::vector<double> x  = linspace(0.0, dVec.maxCoeff());
  std::vector<double> y1 = transform(x, [&](auto x) { return analyticalLoadDisplacementCurve(x); });
  auto p                 = plot(x, y1, dVec, lambdaVec);
  p[0]->line_width(2);
  p[1]->line_width(2);
  p[1]->marker(line_spec::marker_style::asterisk);
  // f->draw();
  //  using namespace std::chrono_literals;
  //  std::this_thread::sleep_for(5s);
}
