#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <vector>
#include <assert.h>

#include "particle.h"
#include "sph_math.h"
#include "eos_base.h"

template <class T>
class ParticleSystem{

private:

    std::vector<Particle<T>> sph_particles;
    int nparticles;
    double h;
    double lambda;
    double nu;
    EOSBase* eos;
    

public:
    ParticleSystem(std::vector<T> r, std::vector<T> v, 
                   std::vector<double> mass, double lh, double lambda, double nu,
                    EOSBase* leos);
    double get_particle_density(int ipart); //< Density around a particle
    double get_density(T ri); //< Density around a point
    void update_acceleration();
    int get_nparticles();
    Particle<T>* get_particle(int ipart);
    //void update_neighbour_table();
    //std::vector<std::vector<int>> get_neighbour_table();
    //std::vector<int> get_neighbour_table(int ipart);
    
};

template <class T>
inline int ParticleSystem<T>::get_nparticles(){return nparticles;};

template <class T>
inline Particle<T>* ParticleSystem<T>::get_particle(int ipart){return &(sph_particles[ipart]); };

template class ParticleSystem<Vec3>;
template class ParticleSystem<Vec2>;

#endif