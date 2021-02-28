#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <vector>
#include <assert.h>

#include "vec3.h"
#include "particle.h"
#include "sph_math.h"
#include "eos_base.h"

class ParticleSystem{

private:

    std::vector<Particle> sph_particles;
    int nparticles;
    double h;
    double lambda;
    double nu;
    EOSBase* eos;

    std::vector<std::vector<int>> neighbour_table;
    

public:
    ParticleSystem(std::vector<Vec3<double>> r, std::vector<Vec3<double>> v, 
                   std::vector<double> mass, double lh,double lambda, double nu,
                    EOSBase* leos);
    double get_particle_density(int ipart); //< Density around a particle
    double get_density(Vec3<double> ri); //< Density around a point
    void update_acceleration();
    int get_nparticles();
    Particle* get_particle(int ipart);
    void update_neighbour_table();
    std::vector<std::vector<int>> get_neighbour_table();
    std::vector<int> get_neighbour_table(int ipart);
    
};

inline int ParticleSystem::get_nparticles(){return nparticles;};
inline Particle* ParticleSystem::get_particle(int ipart){return &(sph_particles[ipart]); };
inline std::vector<std::vector<int>> ParticleSystem::get_neighbour_table(){return neighbour_table;};
#endif