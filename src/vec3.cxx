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

Vec3::Vec3() : Vec2(){ this->z = .0; }

////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor
///
/// Initializes datamembers with the objects x, y and z
/// \param[in] x value desired to the x component
/// \param[in] y value desired to the y component
/// \param[in] z value desired to the z component

Vec3::Vec3(double x, double y, double z) : Vec2(x,y){ this->z = z;}

Vec3 Vec3::operator* (double c){
    Vec3 r(0.,.0,.0);
    r.x = c*this->x;
    r.y = c*this->y;
    r.z = c*this->z;

    return r;
}

Vec3 Vec3::operator+ (const Vec3 v){
    Vec3 r(0.,.0,.0);
    r.x = v.x+this->x;
    r.y = v.y+this->y;
    r.z = v.z+this->z;

    return r;
}