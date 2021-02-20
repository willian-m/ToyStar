#include "integrator_base.h"

IntegratorBase::IntegratorBase(ParticleSystem* lprevious_sys,
                               ParticleSystem* lcurrent_sys,
                               ParticleSystem* lnext_sys,
                               double ldt){
    
    previous_sys = lprevious_sys;
    current_sys = lcurrent_sys;
    next_sys = lnext_sys;
    dt = ldt;
    
    nparticles = current_sys->get_nparticles();
    assert(nparticles == previous_sys->get_nparticles() && "Previous system has different number of particles than the current one");
    assert(nparticles == next_sys->get_nparticles() && "Next system has different number of particles than the current one");

}

IntegratorBase::IntegratorBase(ParticleSystem* lcurrent_sys,
                               ParticleSystem* lnext_sys,
                               double ldt){
    
    previous_sys = NULL;
    current_sys = lcurrent_sys;
    next_sys = lnext_sys;
    dt = ldt;

    nparticles = current_sys->get_nparticles();
    assert(nparticles == next_sys->get_nparticles() && "Next system has different number of particles than the current one");
    
}