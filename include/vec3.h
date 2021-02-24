#ifndef VEC3_H
#define VEC3_H

#include <complex>

template<class T>
class Vec3{

public:
    T x; ///< x coordinate of the vector
    T y; ///< y coordinate of the vector
    T z; ///< z coordinate of the vector

    Vec3();
    Vec3(T x, T y, T z);
    Vec3<T> operator* (double x);
    Vec3<T> operator+ (const Vec3<T> v);
};

//Allowed values of the template
template class Vec3<int>;
template class Vec3<float>;
template class Vec3<double>;
template class Vec3<std::complex<double>>;


#endif