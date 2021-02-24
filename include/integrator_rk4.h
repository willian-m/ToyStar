#ifndef INTEGRATORRK4_H
#define INTEGRATORRK4_H

#include <array>

#include "vec3.h"
#include "particle.h"
#include "integrator_base.h"

class IntegratorRK4 : public IntegratorBase{

private:

    ParticleSystem* buffer;

public:

    IntegratorRK4(ParticleSystem* current, ParticleSystem* next, ParticleSystem* buffer, double dt);
    void do_step();
    void update_system();

};

#endif