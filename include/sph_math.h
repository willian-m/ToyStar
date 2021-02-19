#ifndef SPHMATH_H
#define SPHMATH_H

#include<math.h>

#include "vec3.h"

class SPHMath{
private:
    SPHMath();
    const static double norm3D;

public:
    
    static double kernel_spline3D(Vec3<double> ri, Vec3<double> rj, double h);
    static double distance(Vec3<double> ri, Vec3<double> rj);

};

#endif