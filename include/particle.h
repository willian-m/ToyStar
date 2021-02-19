#ifndef PARTICLE_H
#define PARTICLE_H

#include <complex>
#include "vec3.h"

class Particle{

private:
    Vec3<double> r;      ///< Particle position
    Vec3<double> v;      ///< Particle velocity
    Vec3<double> a;      ///< Particle acceleration
    double mass;         ///< Particle mass
    double pressure;     ///< Pressure in the particle


public:

    //Creates a particle
    Particle(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az, double mass);
    Particle(Vec3<double> r_in, Vec3<double> v_in, Vec3<double> a_in, double mass);

    //Setters

    void set_position(double x, double y, double z);
    void set_velocity(double vx, double vy, double vz);
    void set_acceleration(double ax, double ay, double az);

    void set_position(Vec3<double> x_in);
    void set_velocity(Vec3<double> v_in);
    void set_acceleration(Vec3<double> a_in);
    
    void set_x(double x);
    void set_y(double y);
    void set_z(double z);
    void set_vx(double vx);
    void set_vy(double vy);
    void set_vz(double vz);
    void set_ax(double ax);
    void set_ay(double ay);
    void set_az(double az);

    //Getters
    Vec3<double> get_position();
    Vec3<double> get_velocity();
    Vec3<double> get_acceleration();

    double get_x();
    double get_y();
    double get_z();
    double get_vx();
    double get_vy();
    double get_vz();
    double get_ax();
    double get_ay();
    double get_az();
    double get_mass();
    
};


inline void Particle::set_x(double x){r.x = x;};
inline void Particle::set_y(double y){r.y = y;};
inline void Particle::set_z(double z){r.z = z;};
inline void Particle::set_vx(double vx){v.x = vx;};
inline void Particle::set_vy(double vy){v.y = vy;};
inline void Particle::set_vz(double vz){v.z = vz;};
inline void Particle::set_ax(double ax){a.x = ax;};
inline void Particle::set_ay(double ay){a.y = ay;};
inline void Particle::set_az(double az){a.z = az;};


inline void Particle::set_position(Vec3<double> r_in){this->r = r_in;};
inline void Particle::set_velocity(Vec3<double> v_in){this->v = v_in;};
inline void Particle::set_acceleration(Vec3<double> a_in){this->a = a_in;};



inline Vec3<double> Particle::get_position(){ return r;};
inline Vec3<double> Particle::get_velocity(){ return v;};
inline Vec3<double> Particle::get_acceleration(){ return a;};


inline double Particle::get_x(){ return r.x;};
inline double Particle::get_y(){ return r.y;};
inline double Particle::get_z(){ return r.z;};

inline double Particle::get_vx(){ return v.x;};
inline double Particle::get_vy(){ return v.y;};
inline double Particle::get_vz(){ return v.z;};

inline double Particle::get_ax(){ return a.x;};
inline double Particle::get_ay(){ return a.y;};
inline double Particle::get_az(){ return a.z;};

inline double Particle::get_mass(){ return mass;};

#endif