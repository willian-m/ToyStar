#include "particle_system.h"

template <class T>
ParticleSystem<T>::ParticleSystem(std::vector<T> r, 
                               std::vector<T> v,
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
        sph_particles.push_back( Particle<T>(r[ipart], v[ipart],
                                          T(), mass[ipart]) );
    
    update_acceleration();
}

template <class T>
double ParticleSystem<T>::get_particle_density(int ipart){
    double density = 0;
    T ri = sph_particles[ipart].get_position();
    for (int jpart=0;jpart<nparticles;++jpart){
        density += sph_particles[jpart].get_mass()
                   *SPHMath::kernel_spline(ri,sph_particles[jpart].get_position(),h);
    }
    return density;
}

template <class T>
void ParticleSystem<T>::update_acceleration(){

    //update_neighbour_table();
    for (int ipart=0; ipart<nparticles;++ipart){
        T ri = sph_particles[ipart].get_position();
        T vi = sph_particles[ipart].get_velocity();
        T acc = (-1.)*lambda*ri+(-1.)*nu*vi;
        sph_particles[ipart].set_acceleration(acc);
    }

    for (int ipart=0; ipart<nparticles;++ipart){
        
        T ri = sph_particles[ipart].get_position();
        T acc_i = sph_particles[ipart].get_acceleration();

        double rho_i = get_density(ri);
        double P_i = eos->get_pressure(rho_i);
        double p_over_rho2_i = P_i/(rho_i*rho_i);

        for (int jpart=ipart+1; jpart<nparticles; ++jpart ){
            T rj = sph_particles[jpart].get_position();
            T acc_j = sph_particles[jpart].get_acceleration();
            double m = sph_particles[jpart].get_mass();
            double rho_j = get_density(rj);
            double P_j = eos->get_pressure(rho_j);
            double p_over_rho2_j = P_j/(rho_j*rho_j);
            T gradient = SPHMath::gradient_kernel_spline(ri,rj,h);
            acc_i = acc_i +(-1.)*m*(p_over_rho2_j+p_over_rho2_i)*gradient;
            acc_j = acc_j + m*(p_over_rho2_j+p_over_rho2_i)*gradient;
            sph_particles[jpart].set_acceleration(acc_j);
            sph_particles[ipart].set_acceleration(acc_i);
        }
          
    }
}

template <class T>
double ParticleSystem<T>::get_density(T ri){
    double density = 0;
    for (int jpart =0; jpart<nparticles;++jpart){
        density += sph_particles[jpart].get_mass()
                   *SPHMath::kernel_spline(ri,sph_particles[jpart].get_position(),h);
    }
    return density;
}