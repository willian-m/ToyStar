 // Author: Willian M. Serenone   18/02/21
#include "vec1.h"

/** \class Vec1
 Auxiliary class for three dimensional vectors.
 To create a double vector, use: `Vec3<double> v(1.,2.,3.);`.
 Components can be accessed and modified directly by its datamembers.
 For instance, to update the y-component, use v.y=10;
*/

////////////////////////////////////////////////////////////////////////////////
/// \brief Default constructor
///
/// Initializes datamembers to the default constructor of the template class

Vec1::Vec1(){
    this->x = .0;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor
///
/// Initializes datamembers with the objects x, y and z
/// \param[in] x value desired to the x component
/// \param[in] y value desired to the y component
/// \param[in] z value desired to the z component

Vec1::Vec1(double x){
    this->x = x;
}

Vec1 Vec1::operator*(double c){
    Vec1 r(0.);
    r.x = c*this->x;

    return r;
}

Vec1 Vec1::operator+(const Vec1 v){
    Vec1 r(0.);
    r.x = v.x+this->x;    
    return r;
}