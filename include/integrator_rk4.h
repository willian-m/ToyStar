#ifndef INTEGRATORRK4_H
#define INTEGRATORRK4_H

#include <array>

#include "vec3.h"
#include "particle.h"
#include "integrator_base.h"

class IntegratorRK4 : public IntegratorBase{

private:

    ParticleSystem* buffer;
    void update_system();
    
public:

    IntegratorRK4(ParticleSystem* current, ParticleSystem* next, ParticleSystem* buffer, double dt);
    void do_step();
    

};

#endif