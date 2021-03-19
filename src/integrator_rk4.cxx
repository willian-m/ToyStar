#include "integrator_rk4.h"

template <class T>
void IntegratorRK4<T>::do_step(){
    const double dt_over_two = this->dt*.5;

    T pos;
    T vel;

    std::vector<T> k1_acc, k2_acc, k3_acc, k4_acc;
    std::vector<T> k1_vel, k2_vel, k3_vel, k4_vel;

    //First approximation
    for (int ipart=0; ipart < this->nparticles; ++ipart){
        Particle<T>* particle_current = this->current_sys->get_particle(ipart);
        Particle<T>* particle_buffer = buffer->get_particle(ipart);
        
        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        k1_acc.push_back(particle_current->get_acceleration());
        k1_vel.push_back(vel);
        
        particle_buffer->set_velocity(vel+k1_acc[ipart]*dt_over_two);
        particle_buffer->set_position(pos+k1_vel[ipart]*dt_over_two);
    }
    buffer->update_acceleration(); //Estimation 1 for the system @ t+ dt/2

    //Second approximation
    for (int ipart=0; ipart < this->nparticles; ++ipart){
        Particle<T>* particle_current = this->current_sys->get_particle(ipart);
        Particle<T>* particle_next = this->next_sys->get_particle(ipart);
        Particle<T>* particle_buffer = buffer->get_particle(ipart);
        
        k2_acc.push_back(particle_buffer->get_acceleration());
        k2_vel.push_back(particle_buffer->get_velocity());

        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        particle_next->set_velocity(vel+k2_acc[ipart]*dt_over_two);
        particle_next->set_position(pos+k2_vel[ipart]*dt_over_two);
    }
    this->next_sys->update_acceleration(); //Estimation 2 for the system @ t+ dt/2

    //Third approximation
    for (int ipart=0; ipart < this->nparticles; ++ipart){
        Particle<T>* particle_current = this->current_sys->get_particle(ipart);
        Particle<T>* particle_next = this->next_sys->get_particle(ipart);
        Particle<T>* particle_buffer = buffer->get_particle(ipart);
        
        //First approximation
        k3_acc.push_back(particle_next->get_acceleration());
        k3_vel.push_back(particle_next->get_velocity());
        
        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        particle_buffer->set_velocity(vel+k3_acc[ipart]*this->dt);
        particle_buffer->set_position(pos+k3_vel[ipart]*this->dt);
    }
    buffer->update_acceleration(); //Estimation 1 for the system @ t+ dt

    //Fourth approximation
    for (int ipart=0; ipart < this->nparticles; ++ipart){
        Particle<T>* particle_current = this->current_sys->get_particle(ipart);
        Particle<T>* particle_next = this->next_sys->get_particle(ipart);
        Particle<T>* particle_buffer = buffer->get_particle(ipart);
        
        //First approximation
        k4_acc.push_back(particle_buffer->get_acceleration());
        k4_vel.push_back(particle_buffer->get_velocity());
        
        vel = particle_current->get_velocity();
        pos = particle_current->get_position();

        T avrg_acc = (k1_acc[ipart]+k2_acc[ipart]*2.+k3_acc[ipart]*2.+k4_acc[ipart])*(1./6.);
        T avrg_vel = (k1_vel[ipart]+k2_vel[ipart]*2.+k3_vel[ipart]*2.+k4_vel[ipart])*(1./6.);

        particle_next->set_velocity(vel+avrg_acc*this->dt);
        particle_next->set_position(pos+avrg_vel*this->dt);
        if (SPHMath::distance(particle_current->get_position(),particle_next->get_position()) > this->current_sys->get_h())
            std::cout << "WARNING: Particle displacement bigger than smoothing length." << std::endl;
    }
    this->next_sys->update_acceleration(); //Estimation 2 for the system @ t+ dt
    update_system();

}

template <class T>
void IntegratorRK4<T>::update_system(){

    ParticleSystem<T>* aux = this->current_sys;
    this->current_sys = this->next_sys;
    this->next_sys = aux;
}

template <class T>
IntegratorRK4<T>::IntegratorRK4(ParticleSystem<T>* lcurrent, 
                                ParticleSystem<T>* lnext,
                                ParticleSystem<T>* lbuffer, 
                                double ldt){

    this->current_sys = lcurrent;
    this->next_sys = lnext;
    buffer = lbuffer;
    this->dt = ldt;
    
    this->nparticles = this->current_sys->get_nparticles();
    assert(this->nparticles == buffer->get_nparticles() && "Buffer system has different number of particles than the current one");
    assert(this->nparticles == this->next_sys->get_nparticles() && "Next system has different number of particles than the current one");
};
