#ifndef PARTICLE_H
#define PARTICLE_H

#include <complex>
#include "vec3.h"

template <class T>
class Particle{

private:
    T r;            ///< Particle position
    T v;            ///< Particle velocity
    T a;            ///< Particle acceleration
    double mass;    ///< Particle mass

public:

    //Creates a particle
    //Particle(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az, double mass);
    Particle(T r_in, T v_in, T a_in, double mass);

    //Setters

    //void set_position(double x, double y, double z);
    //void set_velocity(double vx, double vy, double vz);
    //void set_acceleration(double ax, double ay, double az);

    void set_position(T x_in);
    void set_velocity(T v_in);
    void set_acceleration(T a_in);
    void set_mass(double m);
    
    //void set_x(double x);
    //void set_y(double y);
    //void set_z(double z);
    //void set_vx(double vx);
    //void set_vy(double vy);
    //void set_vz(double vz);
    //void set_ax(double ax);
    //void set_ay(double ay);
    //void set_az(double az);

    //Getters
    T get_position();
    T get_velocity();
    T get_acceleration();

    //double get_x();
    //double get_y();
    //double get_z();
    //double get_vx();
    //double get_vy();
    //double get_vz();
    //double get_ax();
    //double get_ay();
    //double get_az();
    double get_mass();
    
};


/*inline void Particle::set_x(double x){r.x = x;};
inline void Particle::set_y(double y){r.y = y;};
inline void Particle::set_z(double z){r.z = z;};
inline void Particle::set_vx(double vx){v.x = vx;};
inline void Particle::set_vy(double vy){v.y = vy;};
inline void Particle::set_vz(double vz){v.z = vz;};
inline void Particle::set_ax(double ax){a.x = ax;};
inline void Particle::set_ay(double ay){a.y = ay;};
inline void Particle::set_az(double az){a.z = az;};*/

template <class T>
inline void Particle<T>::set_position(T r_in){this->r = r_in;};
template <class T>
inline void Particle<T>::set_velocity(T v_in){this->v = v_in;};
template <class T>
inline void Particle<T>::set_acceleration(T a_in){this->a = a_in;};



template <class T>
inline T Particle<T>::get_position(){ return r;};
template <class T>
inline T Particle<T>::get_velocity(){ return v;};
template <class T>
inline T Particle<T>::get_acceleration(){ return a;};


/*inline double Particle::get_x(){ return r.x;};
inline double Particle::get_y(){ return r.y;};
inline double Particle::get_z(){ return r.z;};

inline double Particle::get_vx(){ return v.x;};
inline double Particle::get_vy(){ return v.y;};
inline double Particle::get_vz(){ return v.z;};

inline double Particle::get_ax(){ return a.x;};
inline double Particle::get_ay(){ return a.y;};
inline double Particle::get_az(){ return a.z;};*/

template <class T>
inline double Particle<T>::get_mass(){ return mass;};

template class Particle<Vec3>;
template class Particle<Vec2>;

#endif