// Gmsh project created on Tue Feb 07 09:29:56 2023
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1.0, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface(1) = {1};
