 // Author: Willian M. Serenone   18/02/21
#include "vec3.h"

/** \class Vec3
 Auxiliary class for three dimensional vectors.
 To create a double vector, use: `Vec3<double> v(1.,2.,3.);`.
 Components can be accessed and modified directly by its datamembers.
 For instance, to update the y-component, use v.y=10;
*/

////////////////////////////////////////////////////////////////////////////////
/// \brief Default constructor
///
/// Initializes datamembers to the default constructor of the template class

template<class T>
Vec3<T>::Vec3(){
    this->x = T();
    this->y = T();
    this->z = T();
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor
///
/// Initializes datamembers with the objects x, y and z
/// \param[in] x value desired to the x component
/// \param[in] y value desired to the y component
/// \param[in] z value desired to the z component

template<class T>
Vec3<T>::Vec3(T x, T y, T z){
    this->x = x;
    this->y = y;
    this->z = z;
}



