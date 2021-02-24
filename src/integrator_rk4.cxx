#include "integrator_rk4.h"

void IntegratorRK4::do_step(){
    const double dt_over_two = dt*.5;

    Vec3<double> pos(.0,.0,.0);
    Vec3<double> vel(.0,.0,.0);

    std::vector<Vec3<double>> k1_acc, k2_acc, k3_acc, k4_acc;
    std::vector<Vec3<double>> k1_vel, k2_vel, k3_vel, k4_vel;

    //First approximation
    for (int ipart=0; ipart < nparticles; ++ipart){
        Particle* particle_current = current_sys->get_particle(ipart);
        Particle* particle_buffer = buffer->get_particle(ipart);
        
        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        k1_acc.push_back(particle_current->get_acceleration());
        k1_vel.push_back(vel);
        
        particle_buffer->set_velocity(vel+k1_acc[ipart]*dt_over_two);
        particle_buffer->set_position(pos+k1_vel[ipart]*dt_over_two);
    }
    buffer->update_acceleration(); //Estimation 1 for the system @ t+ dt/2

    //Second approximation
    for (int ipart=0; ipart < nparticles; ++ipart){
        Particle* particle_current = current_sys->get_particle(ipart);
        Particle* particle_next = next_sys->get_particle(ipart);
        Particle* particle_buffer = buffer->get_particle(ipart);
        
        k2_acc.push_back(particle_buffer->get_acceleration());
        k2_vel.push_back(particle_buffer->get_velocity());

        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        particle_next->set_velocity(vel+k2_acc[ipart]*dt_over_two);
        particle_next->set_position(pos+k2_vel[ipart]*dt_over_two);
    }
    next_sys->update_acceleration(); //Estimation 2 for the system @ t+ dt/2

    //Third approximation
    for (int ipart=0; ipart < nparticles; ++ipart){
        Particle* particle_current = current_sys->get_particle(ipart);
        Particle* particle_next = next_sys->get_particle(ipart);
        Particle* particle_buffer = buffer->get_particle(ipart);
        
        //First approximation
        k3_acc.push_back(particle_next->get_acceleration());
        k3_vel.push_back(particle_next->get_velocity());
        
        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        particle_buffer->set_velocity(vel+k3_acc[ipart]*dt);
        particle_buffer->set_position(pos+k3_vel[ipart]*dt);
    }
    buffer->update_acceleration(); //Estimation 1 for the system @ t+ dt

    //Fourth approximation
    for (int ipart=0; ipart < nparticles; ++ipart){
        Particle* particle_current = current_sys->get_particle(ipart);
        Particle* particle_next = next_sys->get_particle(ipart);
        Particle* particle_buffer = buffer->get_particle(ipart);
        
        //First approximation
        k4_acc.push_back(particle_buffer->get_acceleration());
        k4_vel.push_back(particle_buffer->get_velocity());
        
        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        Vec3<double> avrg_acc = (k1_acc[ipart]+k2_acc[ipart]*2.+k3_acc[ipart]*2.+k4_acc[ipart])*(1./6.);
        Vec3<double> avrg_vel = (k1_vel[ipart]+k2_vel[ipart]*2.+k3_vel[ipart]*2.+k4_vel[ipart])*(1./6.);

        particle_next->set_velocity(vel+avrg_acc*dt);
        particle_next->set_position(pos+avrg_vel*dt);
    }
    next_sys->update_acceleration(); //Estimation 2 for the system @ t+ dt

}

void IntegratorRK4::update_system(){

    ParticleSystem* aux = current_sys;
    current_sys = next_sys;
    next_sys = aux;
}

IntegratorRK4::IntegratorRK4(ParticleSystem* lcurrent, 
                             ParticleSystem* lnext,
                             ParticleSystem* lbuffer, 
                             double ldt){

    current_sys = lcurrent;
    next_sys = lnext;
    buffer = lbuffer;
    dt = ldt;
    
    nparticles = current_sys->get_nparticles();
    assert(nparticles == buffer->get_nparticles() && "Buffer system has different number of particles than the current one");
    assert(nparticles == next_sys->get_nparticles() && "Next system has different number of particles than the current one");
};
