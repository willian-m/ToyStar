#include "particle.h"

///TODO: Add functions documentations

/** \class Particle
 A class to store SPH particle properties. There is no dynamic in here
*/
template<class T>
Particle<T>::Particle(T x, T y, T z, T vx, T vy, T vz, T ax, T ay, T az, T mass){

    r.x = x;
    r.y = y;
    r.z = z;

    v.x = vx;
    v.y = vy;
    v.z = vz;

    a.x = ax;
    a.y = ay;
    a.z = az;

    this->mass = mass;
    ///TODO: Init pressure
}

template<class T>
Particle<T>::Particle(Vec3<T> r_in, Vec3<T> v_in, Vec3<T> a_in, T mass){
    r = r_in;
    v = v_in;
    a = a_in;
    this->mass = mass;
    ///TODO: Init pressure
}


template<class T>
void Particle<T>::set_position(T x, T y, T z){
    r.x = x;
    r.y = y;
    r.z = z;
}

template<class T>
void Particle<T>::set_velocity(T vx, T vy, T vz){
    v.x = vx;
    v.y = vy;
    v.z = vz;
}

template<class T>
void Particle<T>::set_acceleration(T ax, T ay, T az){
    a.x = ax;
    a.y = ay;
    a.z = az;
}