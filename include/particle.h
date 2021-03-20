#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <list>
#include "vec2.h"
#include "vec3.h"


template <class T>
class Particle{

private:

    struct NeighbourParticle
    {
        Particle<T>* particle_ptr;
        double distance;
    };
    std::list<NeighbourParticle> neighbour_list;   ///< List of pointers to particles that are in the neighborhood of the particle

    T r;                                                    ///< Particle position
    T v;                                                    ///< Particle velocity
    T a;                                                    ///< Particle acceleration
    double mass;                                            ///< Particle mass
    double density;                                         ///< Particle density
    
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
    void set_density(double rho);

    void add_neighbour_particle(Particle<T>* part, double d);

    int get_num_neighbors();
    auto get_begin_neighbor_list();
    auto get_end_neighbor_list();
    void clear_neighbor_list();
    void erase_neighbor(Particle<T>* part);
    
    //Getters
    T get_position();
    T get_velocity();
    T get_acceleration();

    double get_mass();
    double get_density();
    
};

template <class T>
inline void Particle<T>::set_position(T r_in){this->r = r_in;};
template <class T>
inline void Particle<T>::set_velocity(T v_in){this->v = v_in;};
template <class T>
inline void Particle<T>::set_acceleration(T a_in){this->a = a_in;};
template <class T>
inline void Particle<T>::set_density(double rho){this->density = rho;};



template <class T>
inline T Particle<T>::get_position(){ return r;};
template <class T>
inline T Particle<T>::get_velocity(){ return v;};
template <class T>
inline T Particle<T>::get_acceleration(){ return a;};

template <class T>
inline int Particle<T>::get_num_neighbors(){ return neighbour_list.size();};
template <class T>
inline auto Particle<T>::get_begin_neighbor_list(){ return neighbour_list.begin();};
template <class T>
inline auto Particle<T>::get_end_neighbor_list(){ return neighbour_list.end();};

template <class T>
inline double Particle<T>::get_mass(){ return mass;};
template <class T>
inline double Particle<T>::get_density(){ return density;};

template class Particle<Vec3>;
template class Particle<Vec2>;

#endif