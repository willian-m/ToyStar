#ifndef SPHMATH_H
#define SPHMATH_H

#include<math.h>

#include "vec3.h"

#define TOL 1.E-15

class SPHMath{
private:
    SPHMath();
    const static double norm3D;

public:
    
    static double kernel_spline3D(Vec3<double> ri, Vec3<double> rj, double h);
    static double kernel_spline3D(double distance, double h); //TODO: Implement test
    static double distance(Vec3<double> ri, Vec3<double> rj);
    static Vec3<double> gradient_kernel_spline3D(Vec3<double> ri, Vec3<double> rj, double h); //TODO: Implement test

};

#endif