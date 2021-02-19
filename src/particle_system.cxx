#include "particle_system.h"

ParticleSystem::ParticleSystem(std::vector<Vec3<double>> r, 
                               std::vector<Vec3<double>> v,
                               std::vector<double> mass, double lh
                               ){

    assert(lh > 0);
    h = lh;

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
    std::vector<int> neighbour_table;

    //Make a list of particles in the neighborhood
    for (int ipart=0; ipart<nparticles;++ipart){
        if ( SPHMath::distance(r,sph_particles[ipart].get_position()) < 2*h )
            neighbour_table.push_back(ipart);
    }
    
    for (int ipart : neighbour_table){
        assert(ipart < nparticles);
        density += sph_particles[ipart].get_mass()
                   *SPHMath::kernel_spline3D(r,sph_particles[ipart].get_position(),h);
    }
    return density;
}