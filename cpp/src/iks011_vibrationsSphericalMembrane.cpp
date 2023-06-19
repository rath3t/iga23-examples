// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H

#  include "config.h"

#endif

#include <iostream>
#include <numbers>

#include <dune/grid/uggrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <spdlog/spdlog.h>

#include <ikarus/utils/init.hh>

int main(int argc, char **argv) {
    Ikarus::init(argc, argv);
    using namespace Dune;
    constexpr int gridDim = 2;
    constexpr int dimWorld = 3;
    using Grid = Dune::ALUGrid<gridDim, dimWorld, Dune::simplex, Dune::conforming>;
    auto grid = Dune::GmshReader<Grid>::read("auxiliaryFiles/sphere.msh", false, false);
    grid->globalRefine(0);
    auto gridView = grid->leafGridView();

    spdlog::info("Number of elements: {}", gridView.size(0));
    spdlog::info("Number of vertices: {}", gridView.size(2));

    using namespace Dune::Functions::BasisFactory;
    auto globalbasis = makeBasis(gridView, power<3>(lagrange<1>(), FlatInterleaved()));
    spdlog::info("Number of degrees of freedom: {}", globalbasis.size());

    std::vector<double> values;

    int timeSteps = 200;
    for (int ts = 0; ts < timeSteps; ++ts) {
        auto f = [&](auto x) -> Dune::FieldVector<double, 3> {
            // https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/06%3A_Problems_in_Higher_Dimensions/6.06%3A_Spherically_Symmetric_Vibrations
            double t = (static_cast<double>(ts) / (timeSteps - 1)) * 5.0;
            auto radius = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            auto theta = std::acos(x[2] / radius);
            auto phi = std::copysign(1.0, theta) * std::acos(x[0] / (std::sqrt(x[0] * x[0] + x[1] * x[1])));
            using std::cos;
            using std::sin;
            auto u = 1.0219854770 * sin(3.4641016160 * t) * cos(theta) * sin(theta) * sin(theta) *
                     cos(2.0000000000 * phi);
            Dune::FieldVector<double, 3> UX;
            UX[0] = u * std::sin(theta) * std::cos(phi);
            UX[1] = u * std::sin(theta) * std::sin(phi);
            UX[2] = u * std::cos(theta);
            return UX;
        };

        values.resize(globalbasis.size());
        Dune::Functions::interpolate(globalbasis, values, f);

        /// Write *.vtu files
        auto fh = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(globalbasis, values);

        /// Usage of SubSamplingVTKWriter it subsamples each element and writes this dummy grid to vtk
        Dune::SubsamplingVTKWriter vtkWriter(gridView, Dune::refinementLevels(2), Dune::VTK::conforming);
        vtkWriter.addVertexData(fh, VTK::FieldInfo("discreteVibration", VTK::FieldInfo::Type::vector, 3));
        vtkWriter.write("vibrationSphereInterpolation" + std::to_string(ts));
    }
}
