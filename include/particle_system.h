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

    std::vector<int> get_neighbour_table(Vec3<double> r);
    


public:
    ParticleSystem(std::vector<Vec3<double>> r, std::vector<Vec3<double>> v, 
                   std::vector<double> mass, double lh,double lambda, double nu,
                    EOSBase* leos); //TODO: Test EOS
    double get_density(Vec3<double> r);
    Vec3<double> get_acceleration(int ipart);
    int get_nparticles();
    Particle* get_particle(int ipart);
};

inline int ParticleSystem::get_nparticles(){return nparticles;};
inline Particle* ParticleSystem::get_particle(int ipart){return &(sph_particles[ipart]); };
#endif