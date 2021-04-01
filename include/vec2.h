#ifndef VEC2_H
#define VEC2_H

#include "vec1.h"

class Vec2 : public Vec1{

public:
    double x; ///< x coordinate of the vector
    double y; ///< y coordinate of the vector
    static constexpr double z = 0; ///< We are always constrained to z = 0 plane

    Vec2();
    Vec2(double x, double y);
    Vec2 operator* (double x);
    Vec2 operator+ (const Vec2 v);
};

inline Vec2 operator*(double x, Vec2 v){ return v*x; };

#endif