#include "particle_system.h"

ParticleSystem::ParticleSystem(std::vector<Vec3<double>> r, 
                               std::vector<Vec3<double>> v,
                               std::vector<double> mass, double lh,
                               double llambda, double lnu,
                               EOSBase* leos
                               ){

    assert(lh > 0);
    h = lh;
    eos = leos;
    assert(lnu > 0);
    assert(llambda > 0);

    nu = lnu;
    lambda = llambda;

    double n_r = r.size();
    int n_v = v.size();
    int n_m = mass.size();

    assert(n_r == n_v && "ERROR: Position and velocity arrays have different sizes");
    assert(n_v == n_m && "ERROR: Velocity and mass arrays have different sizes");

    


    nparticles = n_r;

    for(int ipart=0; ipart< nparticles; ++ipart)
        sph_particles.push_back( Particle(r[ipart], v[ipart],
                                          Vec3<double>(.0,.0,.0), mass[ipart]) );
    

}

double ParticleSystem::get_density(Vec3<double> r){
    double density = 0;
    
    std::vector<int> neighbour_table = get_neighbour_table(r);
    
    for (int ipart : neighbour_table){
        density += sph_particles[ipart].get_mass()
                   *SPHMath::kernel_spline3D(r,sph_particles[ipart].get_position(),h);
    }
    return density;
}

Vec3<double> ParticleSystem::get_acceleration(int ipart){

    Vec3<double> ri = sph_particles[ipart].get_position();
    std::vector<int> neighbour_table = get_neighbour_table(ri);

    double rho_i = get_density(ri);
    double P_i = eos->get_pressure(rho_i);
    double p_over_rho2_i = P_i/(rho_i*rho_i);

    Vec3<double> acc = ri*(-1.*lambda)+sph_particles[ipart].get_velocity()*(-1.*nu);
    
    
    for (int jpart : neighbour_table){
        double m = sph_particles[jpart].get_mass();
        Vec3<double> r_j = sph_particles[jpart].get_position();
        double rho_j = get_density(r_j);
        double P_j = eos->get_pressure(rho_i);
        double p_over_rho2_j = P_i/(rho_j*rho_j);
        Vec3<double> gradient = SPHMath::gradient_kernel_spline3D(ri,r_j,this->h);
        acc = acc + gradient*p_over_rho2_j;   
    }

    return acc;
}

std::vector<int> ParticleSystem::get_neighbour_table(Vec3<double> r){
    std::vector<int> neighbour_table;
    for (int ipart=0; ipart<nparticles;++ipart){
        if ( SPHMath::distance(r,sph_particles[ipart].get_position()) < 2*h )
        neighbour_table.push_back(ipart);
    }
    return neighbour_table;
}