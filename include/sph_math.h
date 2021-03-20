#ifndef SPHMATH_H
#define SPHMATH_H

#include<math.h>
#include<typeinfo>

#include "vec3.h"

#define TOL 1.E-15

class SPHMath{
private:
    SPHMath();
    const static double norm3D;
    const static double norm2D;

public:
    
    static double kernel_spline(Vec3 ri, Vec3 rj, double h);
    static double kernel_spline(Vec2 ri, Vec2 rj, double h);
    static double kernel_spline3D(double distance, double h);
    static double kernel_spline2D(double distance, double h);
    static double distance(Vec3 ri, Vec3 rj);
    static double distance(Vec2 ri, Vec2 rj);
    static Vec3 gradient_kernel_spline(Vec3 ri, Vec3 rj, double h);
    static Vec2 gradient_kernel_spline(Vec2 ri, Vec2 rj, double h);
    template <class T>
    static double kernel_spline(double d, double h);

};

#endif