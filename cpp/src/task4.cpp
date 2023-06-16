// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/common/parametertreeparser.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/enhancedAssumedStrains.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/duneUtilities.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

//using namespace Ikarus;
//using namespace Dune::Indices;

int main(int argc, char **argv) {
    Ikarus::init(argc, argv);

    constexpr int gridDim = 2;
    double lambdaLoad = 1;
    constexpr int basis_order = 1;

    /// Read input parameters

    Dune::ParameterTree parameterSet;
    Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

    const Dune::ParameterTree &gridParameters = parameterSet.sub("GridParameters");
    const Dune::ParameterTree &controlParameters = parameterSet.sub("ControlParameters");
    const Dune::ParameterTree &materialParameters = parameterSet.sub("MaterialParameters");
    const Dune::ParameterTree &elementParameters = parameterSet.sub("ElementParameters");

    const auto E = materialParameters.get<double>("E");
    const auto nu = materialParameters.get<double>("nu");
    const auto refinement_level = gridParameters.get<int>("refinement");
    const auto boundaryLoadPredicate = gridParameters.get<std::string>("boundaryloadpredicate");
    const auto gridFileName = gridParameters.get<std::string>("filename");
    const auto numberOfEASParameters = elementParameters.get<int>("numberOfEASParameters");

    ///TODO: Read grid file from folder auxiliaryFiles/ with the name gridFileName and create a UGGrid


    grid->globalRefine(refinement_level);
    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    ///TODO: Create a lagrange basis for the 2D case
    ///see Dune page 336 Function space bases

    ///TODO: Clamp the left boundary

    ///TODO: Declare a vector "fes" of linear elastic 2D planar solid elements decorated with EAS

    ///TODO: Function for volume load- here: returns zero


    ///TODO: Neumann boundary load in vertical direction


    /// Python function which could be used to obtain the vertices at the right edge
    std::string lambdaNeumannVertices = std::string(boundaryLoadPredicate);
    Python::start();
    Python::Module main = Python::import("__main__");
    Python::run("import math");

    Python::runStream() << std::endl << "import sys" << std::endl << "import os" << std::endl;

    const auto &indexSet = gridView.indexSet();

    /// Flagging the vertices on which neumann load is applied as true
    Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
    auto pythonNeumannVertices = Python::make_function<bool>(Python::evaluate(lambdaNeumannVertices));
    for (auto &&vertex: vertices(gridView)) {
        bool isNeumann = pythonNeumannVertices(vertex.geometry().corner(0));
        neumannVertices[indexSet.index(vertex)] = isNeumann;
    }

    BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

    ///TODO: Add the linear elastic 2D planar solid elements decorated with EAS to the vector "fes"
    for (auto &element: elements(gridView)) {
        auto localView = basis->localView();


    }

    ///TODO: Create a sparse assembler

    ///TODO: Define "elastoStatics" affordances and create functions for stiffness matrix and residual calculations


    Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis->size());

    ///TODO: Create a non-linear operator
    auto startAssembly = std::chrono::high_resolution_clock::now();



    auto stopAssembly = std::chrono::high_resolution_clock::now();
    auto durationAssembly = duration_cast<std::chrono::milliseconds>(stopAssembly - startAssembly);
    spdlog::info("The assembly took {:>6d} milliseconds with {:>7d} dofs",
                 durationAssembly.count(), basis->size());

    const auto &K = nonLinOp.derivative();
    const auto &Fext = nonLinOp.value();

    ///TODO: solve the linear system

    auto startSolver = std::chrono::high_resolution_clock::now();


    auto stopSolver = std::chrono::high_resolution_clock::now();
    auto durationSolver = duration_cast<std::chrono::milliseconds>(stopSolver - startSolver);
    spdlog::info("The solver took {} milliseconds ",
                 durationSolver.count());

    ///TODO: Postprocess


    Dune::VTKWriter vtkWriter(gridView, Dune::VTK::conforming);
    vtkWriter.addVertexData(dispGlobalFunc,
                            Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
    vtkWriter.write(gridFileName);
    std::cout << std::endl << "################################################################" << std::endl;
}
