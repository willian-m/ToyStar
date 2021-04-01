#ifndef INTEGRATORRK4_H
#define INTEGRATORRK4_H

#include <array>
#include <iostream>

#include "vec1.h"
#include "vec2.h"
#include "vec3.h"
#include "particle.h"
#include "integrator_base.h"

template <class T>
class IntegratorRK4 : public IntegratorBase<T>{

private:

    ParticleSystem<T>* buffer;
    void update_system();
    
public:

    IntegratorRK4(ParticleSystem<T>* current, 
                  ParticleSystem<T>* next,
                  ParticleSystem<T>* buffer, double dt);
    void do_step();
    

};

template class IntegratorRK4<Vec3>;
template class IntegratorRK4<Vec2>;
template class IntegratorRK4<Vec1>;

#endif