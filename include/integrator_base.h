//Abstract class defining how a common interface for all integrators 

#ifndef IntegratorBase_H
#define IntegratorBase_H

#include "particle_system.h"

template <class T>
class IntegratorBase{

protected:
    ParticleSystem<T>* previous_sys;
    ParticleSystem<T>* current_sys;
    ParticleSystem<T>* next_sys;
    double dt;
    int nparticles;
    virtual void update_system() = 0;

public:
    
    virtual void do_step() = 0; ///< Derived class must specify how to perform time-steps
};

template class IntegratorBase<Vec3>;
template class IntegratorBase<Vec2>;
template class IntegratorBase<Vec1>;

#endif