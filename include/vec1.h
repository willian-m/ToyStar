#ifndef VEC1_H
#define VEC1_H

class Vec1{

public:
    double x; ///< x coordinate of the vector
    static constexpr double y = 0; ///<  We are always constrained to y = 0
    static constexpr double z = 0; ///< We are always constrained to z = 0

    Vec1();
    Vec1(double x);
    Vec1 operator* (double x);
    Vec1 operator+ (const Vec1 v);
};

inline Vec1 operator*(double x, Vec1 v){ return v*x; };

#endif