#ifndef VEC2_H
#define VEC2_H

#include <complex>

class Vec2{

public:
    double x; ///< x coordinate of the vector
    double y; ///< y coordinate of the vector

    Vec2();
    Vec2(double x, double y);
    Vec2 operator* (double x);
    Vec2 operator+ (const Vec2 v);
};

inline Vec2 operator*(double x, Vec2 v){ return v*x; };

#endif