// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H

#  include "config.h"

#endif

#include <iostream>
#include <numbers>

#include <dune/alugrid/grid.hh>
#include <dune/common/exceptions.hh>          // We use exceptions
#include <dune/common/parallel/mpihelper.hh>  // An initializer of MPI
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <dune/common/parametertreeparser.hh>

#define GridDim 3

int main(int argc, char **argv) {
    Dune::MPIHelper::instance(argc, argv);
    Dune::ParameterTree parameterSet;
    Dune::ParameterTreeParser::readINITree("/workspaces/iga23-examples/cpp/src/auxiliaryFiles/parameterTask1.parset", parameterSet);

    const int refinementLevels = parameterSet.get<int>("refinementLevels");
    using namespace Dune;

#if GridDim == 2
    Eigen::Vector2d v(1, 0);
    std::array<FieldVector<double, 2>, 6> corners0;
    Eigen::Rotation2D<double> R;
    R.angle() = 0.0;
    for (auto &corner: corners0) {
      Eigen::Vector2d a = R * v;
      corner[0] = a[0];
      corner[1] = a[1];
      R.angle() += 60.0 / 180.0 * std::numbers::pi;
    }

    constexpr int gridDim = 2;
    constexpr int dimWorld = 2;

    Dune::GridFactory<Dune::ALUGrid<gridDim, dimWorld, Dune::simplex, Dune::conforming>> gridFactory;
    gridFactory.insertVertex({0, 0});
    gridFactory.insertVertex(corners0[0]);
    gridFactory.insertVertex(corners0[1]);
    gridFactory.insertVertex(corners0[2]);
    gridFactory.insertVertex(corners0[3]);
    gridFactory.insertVertex(corners0[4]);
    gridFactory.insertVertex(corners0[5]);

    gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 1, 2}); // ordering see dune book page 128 Figure 5.12
    gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 2, 3});
    gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 3, 4});
    gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 4, 5});
    gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 5, 6});
    gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 6, 1});

#elif GridDim == 3
    constexpr int gridDim = 3;
    constexpr int dimWorld = 3;
    Dune::GridFactory<Dune::ALUGrid<gridDim, dimWorld, Dune::simplex, Dune::conforming>> gridFactory;

    gridFactory.insertVertex({0, 0, 0});
    gridFactory.insertVertex({1, 0, 0});
    gridFactory.insertVertex({0.5, 1, 0});
    gridFactory.insertVertex({0.5, 0.5, 1});

    gridFactory.insertElement(Dune::GeometryTypes::tetrahedron,
                              {0, 1, 2, 3}); // ordering see dune book page 128 Figure 5.13
#endif

    auto grid = gridFactory.createGrid();
    grid->globalRefine(refinementLevels);
    auto gridView = grid->leafGridView();
    auto &indexSet = gridView.indexSet();

    std::vector<double> elementmarker;
    const int numberOfElements = gridView.size(0);
    elementmarker.resize(numberOfElements);
    for (const auto &ele: elements(gridView)) {
        auto elementId = indexSet.index(ele);
        if (ele.hasBoundaryIntersections())
            elementmarker[elementId] = 1;
        else
            elementmarker[elementId] = 0;
    }

    Dune::VTKWriter vtkWriter(gridView);
    vtkWriter.addCellData(elementmarker, "isAtBoundary");
    vtkWriter.write("HexagonWithFlaggedBoundary");

}
