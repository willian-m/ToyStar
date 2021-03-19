#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "vec2.h"
#include "vec3.h"

template <class T>
class Particle{

private:
    T r;                                        ///< Particle position
    T v;                                        ///< Particle velocity
    T a;                                        ///< Particle acceleration
    double mass;                                ///< Particle mass
    std::vector<Particle<T>*> neighbour_list;   ///< List of pointers to particles that are in the neighborhood of the particle
    std::vector<double> neighbour_list_distance;///< List with distance to neighbour particles
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
    void add_neighbour_particle(Particle<T>* part, double d);
    int get_num_neighbors();
    Particle<T>* get_neighbor(int i);
    double get_neighbor_distance(int i);
    void clear_neighbors();
    
    //Getters
    T get_position();
    T get_velocity();
    T get_acceleration();

    double get_mass();
    
};

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

template <class T>
inline int Particle<T>::get_num_neighbors(){ return neighbour_list.size();};
template <class T>
inline Particle<T>* Particle<T>::get_neighbor(int i){ return neighbour_list[i];};
template <class T>
inline double Particle<T>::get_neighbor_distance(int i){ return neighbour_list_distance[i];};

template <class T>
inline double Particle<T>::get_mass(){ return mass;};

template class Particle<Vec3>;
template class Particle<Vec2>;

#endif