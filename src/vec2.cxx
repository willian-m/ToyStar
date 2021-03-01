 // Author: Willian M. Serenone   18/02/21
#include "vec2.h"

/** \class Vec2
 Auxiliary class for three dimensional vectors.
 To create a double vector, use: `Vec3<double> v(1.,2.,3.);`.
 Components can be accessed and modified directly by its datamembers.
 For instance, to update the y-component, use v.y=10;
*/

////////////////////////////////////////////////////////////////////////////////
/// \brief Default constructor
///
/// Initializes datamembers to the default constructor of the template class

Vec2::Vec2(){
    this->x = .0;
    this->y = .0;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor
///
/// Initializes datamembers with the objects x, y and z
/// \param[in] x value desired to the x component
/// \param[in] y value desired to the y component
/// \param[in] z value desired to the z component

Vec2::Vec2(double x, double y){
    this->x = x;
    this->y = y;
}

Vec2 Vec2::operator*(double c){
    Vec2 r(0.,.0);
    r.x = c*this->x;
    r.y = c*this->y;

    return r;
}

Vec2 Vec2::operator+(const Vec2 v){
    Vec2 r(0.,.0);
    r.x = v.x+this->x;
    r.y = v.y+this->y;
    
    return r;
}