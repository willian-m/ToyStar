//Abstract class defining how a common interface for all integrators 

#ifndef IntegratorBase_H
#define IntegratorBase_H

#include "particle_system.h"

class IntegratorBase{

protected:
    ParticleSystem* previous_sys;
    ParticleSystem* current_sys;
    ParticleSystem* next_sys;
    double dt;
    int nparticles;
    virtual void update_system() = 0;

public:
    /*IntegratorBase(ParticleSystem* previous_sys,
                   ParticleSystem* current_sys,
                   ParticleSystem* next_sys,
                   double dt);

    IntegratorBase(ParticleSystem* current_sys,
                   ParticleSystem* next_sys,
                   double dt);*/

    virtual void do_step() = 0; ///< Derived class must specify how to perform time-steps
    

    
};

#endif