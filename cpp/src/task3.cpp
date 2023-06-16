// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H

#  include "config.h"

#endif

#include <iostream>
#include <numbers>

#include <dune/alugrid/grid.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/utils/init.hh>

int main(int argc, char **argv) {
  Ikarus::init(argc, argv);
  std::cout << "Hello to Ikarus exercises Task 3." << std::endl;
  
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const int refinementLevels = parameterSet.get<int>("refinementLevels");
  using namespace Dune;
  /// 2D grid
  constexpr int gridDim  = 2;
  constexpr int dimWorld = 2;
  Dune::GridFactory<Dune::ALUGrid<gridDim, dimWorld, Dune::simplex, Dune::conforming>> gridFactory;
  Eigen::Vector2d v(1, 0);
  std::array<FieldVector<double, 2>, 6> corners0;
  Eigen::Rotation2D<double> R;
  R.angle() = 0.0;
  for (auto &corner : corners0) {
    Eigen::Vector2d a = R * v;
    corner[0]         = a[0];
    corner[1]         = a[1];
    R.angle() += 60.0 / 180.0 * std::numbers::pi;
  }
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex(corners0[0]);
  gridFactory.insertVertex(corners0[1]);
  gridFactory.insertVertex(corners0[2]);
  gridFactory.insertVertex(corners0[3]);
  gridFactory.insertVertex(corners0[4]);
  gridFactory.insertVertex(corners0[5]);

  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 1, 2});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 2, 3});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 3, 4});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 4, 5});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 5, 6});
  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 6, 1});

  auto grid = gridFactory.createGrid();
  grid->globalRefine(refinementLevels);
  auto gridView  = grid->leafGridView();
  auto &indexSet = gridView.indexSet();
  using GridView = decltype(gridView);

  spdlog::info("Number of elements: {}", gridView.size(0));
  spdlog::info("Number of vertices: {}", gridView.size(2));

  std::vector<double> values;
  double t = 2.9;
  auto f   = [&](auto x) -> double {
    // https://en.wikipedia.org/wiki/Vibrations_of_a_circular_membrane
    return std::sin(8.65373 * t) * std::cyl_bessel_j(0, 8.65373 * std::sqrt(x[0] * x[0] + x[1] * x[1]));
  };

  /// TODO Create a first-order Lagrangian basis on the leaf grid view

  /// TODO Create a vector, where the interpolated w(xi , t) should be stored

  /// TODO Create a local view of the Lagrangian basis

  /// Loop over all elements
  for (const auto &ele : elements(gridView)) {
    /// TODO Bind local view to current grid element

    /// TODO access local Lagrangian basis from local view

    /// TODO get "finite element" from local Lagrangian basis

    /// TODO Create geometry object from gridElement

    /// Loop over vertices of grid elements
    for (int i = 0; i < geo.corners(); ++i) {
      auto coordsCorner = geo.corner(i);

      /// TODO Evaluate w(xi , t)

      /// TODO get globalIndex where w(xi,t) should be stored

      /// TODO store result in global vector
    }

    /// TODO Create a DiscreteGlobalBasisFunction from the global basis and the vector of coefficients

    /// TODO Write solution to vtk
    Dune::VTKWriter vtkWriter(gridView);
  }
}
