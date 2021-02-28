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
    assert(lnu >= 0.);
    assert(llambda >= 0.);

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
    
    update_acceleration();
}

double ParticleSystem::get_particle_density(int ipart){
    double density = 0;
    Vec3<double> ri = sph_particles[ipart].get_position();
    for (int jpart : neighbour_table[ipart]){
        density += sph_particles[jpart].get_mass()
                   *SPHMath::kernel_spline3D(ri,sph_particles[ipart].get_position(),h);
    }
    return density;
}

void ParticleSystem::update_acceleration(){

    //update_neighbour_table();
    for (int ipart=0; ipart<nparticles;++ipart){
        Vec3<double> ri = sph_particles[ipart].get_position();
        Vec3<double> vi = sph_particles[ipart].get_velocity();
        Vec3<double> acc = ri*(-1.*lambda)+vi*(-1.*nu);
        sph_particles[ipart].set_acceleration(acc);
    }

    for (int ipart=0; ipart<nparticles;++ipart){
        
        Vec3<double> ri = sph_particles[ipart].get_position();
        Vec3<double> acc_i = sph_particles[ipart].get_acceleration();

        double rho_i = get_density(ri);
        double P_i = eos->get_pressure(rho_i);
        double p_over_rho2_i = P_i/(rho_i*rho_i);

        for (int jpart=ipart+1; jpart<nparticles; ++jpart ){
            Vec3<double> rj = sph_particles[jpart].get_position();
            Vec3<double> acc_j = sph_particles[jpart].get_acceleration();
            double m = sph_particles[jpart].get_mass();
            double rho_j = get_density(rj);
            double P_j = eos->get_pressure(rho_j);
            double p_over_rho2_j = P_j/(rho_j*rho_j);
            Vec3<double> gradient = SPHMath::gradient_kernel_spline3D(ri,rj,h);
            acc_i = acc_i + gradient*(p_over_rho2_j+p_over_rho2_i)*(-1.*m);
            acc_j = acc_j + gradient*(p_over_rho2_j+p_over_rho2_i)*m;
            sph_particles[jpart].set_acceleration(acc_j);
            sph_particles[ipart].set_acceleration(acc_i);
        }
        
        
    }
}

void ParticleSystem::update_neighbour_table(){
    neighbour_table.clear();
    std::vector<int> neighbors = {};
    for (int ipart=0; ipart<nparticles;++ipart){
        neighbors = get_neighbour_table(ipart);
        neighbour_table.push_back(neighbors);
    }
}

std::vector<int> ParticleSystem::get_neighbour_table(int ipart){
    Vec3<double> ri = sph_particles[ipart].get_position();
    std::vector<int> neighbour = {};
    for (int jpart=0; jpart<nparticles;++jpart) {
        Vec3<double> rj = sph_particles[jpart].get_position();
        if (SPHMath::kernel_spline3D(ri,rj,h) > TOL ) neighbour.push_back(jpart);
    };
    return neighbour;
}
double ParticleSystem::get_density(Vec3<double> ri){
    double density = 0;
    for (int jpart =0; jpart<nparticles;++jpart){
        density += sph_particles[jpart].get_mass()
                   *SPHMath::kernel_spline3D(ri,sph_particles[jpart].get_position(),h);
    }
    return density;
}