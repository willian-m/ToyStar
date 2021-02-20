#include "integrator_rk4.h"

void IntegratorRK4::do_step(){
    const double dt_over_two = dt*.5;
    for (int ipart=0; ipart < nparticles; ++ipart){

        Particle* particle_current = current_sys->get_particle(ipart);
        Particle* particle_next = next_sys->get_particle(ipart);
        
        //First approximation
        std::array<Vec3<double>,2> k1 = get_k1(ipart);
        particle_next->set_velocity(particle_current->get_velocity()+k1[1]*dt_over_two);
        particle_next->set_position(particle_current->get_position()+k1[0]*dt_over_two);

        //Second approximation
        std::array<Vec3<double>,2> k2 = get_k2(ipart);
        particle_next->set_velocity(particle_current->get_velocity()+k2[1]*dt_over_two);
        particle_next->set_position(particle_current->get_position()+k2[0]*dt_over_two);
        
        //Third approximation
        std::array<Vec3<double>,2> k3 = get_k3(ipart);
        particle_next->set_velocity(particle_current->get_velocity()+k3[1]*dt);
        particle_next->set_position(particle_current->get_position()+k3[0]*dt);

        //Final approximation
        std::array<Vec3<double>,2> k4 = get_k4(ipart);
        particle_next->set_velocity(
            particle_current->get_velocity() + (k1[1] + k2[1]*2. + k3[1]*2. + k4[1])*(dt/6.) );
        particle_next->set_position(
            particle_current->get_position() + (k1[0] + k2[0]*2. + k3[0]*2. + k4[0])*(dt/6.) );
    }

}

void IntegratorRK4::update_system(){

    ParticleSystem* aux = current_sys;
    current_sys = next_sys;
    next_sys = aux;
}

std::array<Vec3<double>,2> IntegratorRK4::get_k1(int ipart){
    std::array<Vec3<double>,2> k1;
    k1[0] = current_sys->get_particle(ipart)->get_velocity();
    k1[1] = current_sys->get_acceleration(ipart);
    return k1;
}

std::array<Vec3<double>,2> IntegratorRK4::get_k2(int ipart){
    std::array<Vec3<double>,2> k2;
    k2[0] = next_sys->get_particle(ipart)->get_velocity();
    k2[1] = next_sys->get_acceleration(ipart);
    return k2;
}

std::array<Vec3<double>,2> IntegratorRK4::get_k3(int ipart){
    std::array<Vec3<double>,2> k3;
    k3[0] = next_sys->get_particle(ipart)->get_velocity();
    k3[1] = next_sys->get_acceleration(ipart);
    return k3;
}

std::array<Vec3<double>,2> IntegratorRK4::get_k4(int ipart){
    std::array<Vec3<double>,2> k4;
    k4[0] = next_sys->get_particle(ipart)->get_velocity();
    k4[1] = next_sys->get_acceleration(ipart);
    return k4;
}