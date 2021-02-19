#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <vector>
#include <assert.h>

#include "vec3.h"
#include "particle.h"
#include "sph_math.h"

class ParticleSystem{

private:

    std::vector<Particle> sph_particles;
    int nparticles;
    double h;

public:
    ParticleSystem(std::vector<Vec3<double>> r, std::vector<Vec3<double>> v, 
                   std::vector<double> mass, double lh);
    double get_density(Vec3<double> r);
};

#endif