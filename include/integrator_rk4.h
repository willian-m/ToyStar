#ifndef INTEGRATORRK4_H
#define INTEGRATORRK4_H

#include <array>

#include "vec3.h"
#include "particle.h"
#include "integrator_base.h"

class IntegratorRK4 : public IntegratorBase{

private:

    std::array<Vec3<double>,2> get_k1(int ipart);
    std::array<Vec3<double>,2> get_k2(int ipart);
    std::array<Vec3<double>,2> get_k3(int ipart);
    std::array<Vec3<double>,2> get_k4(int ipart);


public:

    IntegratorRK4(ParticleSystem* current, ParticleSystem* next, double dt);
    void do_step();
    void update_system();

};

inline IntegratorRK4::IntegratorRK4(ParticleSystem* current, 
                                    ParticleSystem* next, 
                                    double dt):IntegratorBase(current, next, dt){};

#endif