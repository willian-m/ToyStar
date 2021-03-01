#include "particle.h"

///TODO: Add functions documentations

/** \class Particle
 A class to store SPH particle properties. There is no dynamic in here
*/

/*Particle::Particle(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az, double mass){

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

}*/

template <class T>
Particle<T>::Particle(T r_in, T v_in, T a_in, double mass){
    r = r_in;
    v = v_in;
    a = a_in;
    this->mass = mass;
}



/*void Particle::set_position(double x, double y, double z){
    r.x = x;
    r.y = y;
    r.z = z;
}


void Particle::set_velocity(double vx, double vy, double vz){
    v.x = vx;
    v.y = vy;
    v.z = vz;
}


void Particle::set_acceleration(double ax, double ay, double az){
    a.x = ax;
    a.y = ay;
    a.z = az;
}*/