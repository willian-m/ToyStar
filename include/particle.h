#ifndef PARTICLE_H
#define PARTICLE_H

#include <complex>
#include "vec3.h"

template<class T>
class Particle{

private:
    Vec3<T> r;      ///< Particle position
    Vec3<T> v;      ///< Particle velocity
    Vec3<T> a;      ///< Particle acceleration
    T mass;         ///< Particle mass
    T pressure;     ///< Pressure in the particle


public:

    //Creates a particle
    Particle(T x, T y, T z, T vx, T vy, T vz, T ax, T ay, T az, T mass);
    Particle(Vec3<T> r_in, Vec3<T> v_in, Vec3<T> a_in, T mass);

    //Setters

    void set_position(T x, T y, T z);
    void set_velocity(T vx, T vy, T vz);
    void set_acceleration(T ax, T ay, T az);

    void set_position(Vec3<T> x_in);
    void set_velocity(Vec3<T> v_in);
    void set_acceleration(Vec3<T> a_in);
    
    void set_x(T x);
    void set_y(T y);
    void set_z(T z);
    void set_vx(T vx);
    void set_vy(T vy);
    void set_vz(T vz);
    void set_ax(T ax);
    void set_ay(T ay);
    void set_az(T az);

    //Getters
    Vec3<T> get_position();
    Vec3<T> get_velocity();
    Vec3<T> get_acceleration();

    T get_x();
    T get_y();
    T get_z();
    T get_vx();
    T get_vy();
    T get_vz();
    T get_ax();
    T get_ay();
    T get_az();
    
};

//Allowed values of the template
template class Particle<int>;
template class Particle<float>;
template class Particle<double>;
template class Particle<std::complex<double>>;

template<class T>
inline void Particle<T>::set_x(T x){r.x = x;};
template<class T>
inline void Particle<T>::set_y(T y){r.y = y;};
template<class T>
inline void Particle<T>::set_z(T z){r.z = z;};
template<class T>
inline void Particle<T>::set_vx(T vx){v.x = vx;};
template<class T>
inline void Particle<T>::set_vy(T vy){v.y = vy;};
template<class T>
inline void Particle<T>::set_vz(T vz){v.z = vz;};
template<class T>
inline void Particle<T>::set_ax(T ax){a.x = ax;};
template<class T>
inline void Particle<T>::set_ay(T ay){a.y = ay;};
template<class T>
inline void Particle<T>::set_az(T az){a.z = az;};

template<class T>
inline void Particle<T>::set_position(Vec3<T> r_in){this->r = r_in;};
template<class T>
inline void Particle<T>::set_velocity(Vec3<T> v_in){this->v = v_in;};
template<class T>
inline void Particle<T>::set_acceleration(Vec3<T> a_in){this->a = a_in;};


template<class T>
inline Vec3<T> Particle<T>::get_position(){ return r;};
template<class T>
inline Vec3<T> Particle<T>::get_velocity(){ return v;};
template<class T>
inline Vec3<T> Particle<T>::get_acceleration(){ return a;};

template<class T>
T Particle<T>::get_x(){ return r.x;};
template<class T>
T Particle<T>::get_y(){ return r.y;};
template<class T>
T Particle<T>::get_z(){ return r.z;};

template<class T>
T Particle<T>::get_vx(){ return v.x;};
template<class T>
T Particle<T>::get_vy(){ return v.y;};
template<class T>
T Particle<T>::get_vz(){ return v.z;};

template<class T>
T Particle<T>::get_ax(){ return a.x;};
template<class T>
T Particle<T>::get_ay(){ return a.y;};
template<class T>
T Particle<T>::get_az(){ return a.z;};

#endif