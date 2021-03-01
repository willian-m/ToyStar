#ifndef VEC3_H
#define VEC3_H

#include "vec2.h"

class Vec3 : public Vec2 {

public:
    double z; ///< z coordinate of the vector

    Vec3();
    Vec3(double x, double y, double z);
    Vec3 operator* (double x);
    Vec3 operator+ (const Vec3 v);
};

inline Vec3 operator* (double x, Vec3 v){ return v*x; };

#endif