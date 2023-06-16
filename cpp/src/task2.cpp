// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <Eigen/Core>
#include <Eigen/Dense>

#define GridDim 2

int main(int argc, char **argv) {
    Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello to Ikarus exercises Task 2." << std::endl;

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

    /// use gridFactory to insert vertices and elements
    

#else

    //TODO: 3D tetrahedron


#endif

    /// create grid from gridFactory


    /// loop over elements and mark those on the boundary


    /// write result to file for Paraview


}
