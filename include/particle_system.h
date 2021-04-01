#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <vector>

#include <assert.h>
#include <algorithm>
#include <typeinfo>

#include "vec1.h"
#include "vec2.h"
#include "vec3.h"
#include "sph_math.h"
#include "eos_base.h"
#include "particle.h"
#include "cell.h"

template <class T>
class ParticleSystem{

private:

    std::vector<Particle<T>> sph_particles;
    int nparticles;//< Number of particles in the system
    double h;      //< SPH smearing parameter  
    double lambda; //< Gravity 
    double nu;     //< Dampening  
    EOSBase* eos;  //< Pointer to EOS to be used

    //Grid variables
    std::vector<Cell<T>> grid;  //< Grid where particles are organized
    int nx, ny, nz; //< Number of cells in each side of the grid
    int ncells;     //< Total number of cells in the grid
    double lim_x[2], lim_y[2], lim_z[2];
    

    //Auxiliary functions
    void build_cell_neighbour_list();
    void add_particles_to_grid();
    int get_grid_idx(Vec3 r);
    int get_grid_idx(Vec2 r);
    int get_grid_idx(Vec1 r);
    void clear_grid();
    void update_densities();


public:
    ParticleSystem(std::vector<T> r, std::vector<T> v, 
                   std::vector<double> mass, double lh, double lambda, double nu,
                    EOSBase* leos);
    double get_point_density(T ri); //< Density around a point
    void update_acceleration();
    int get_nparticles();
    double get_h();
    Particle<T>* get_particle(int ipart);

    
    void setup_grid(Vec3 lower_bound, Vec3 upper_bounds);
    void find_neighbor_particles();
    
};

template <class T>
inline int ParticleSystem<T>::get_nparticles(){return nparticles;};

template <class T>
inline Particle<T>* ParticleSystem<T>::get_particle(int ipart){return &(sph_particles[ipart]); };
template <class T>
inline double ParticleSystem<T>::get_h(){return h;};


template class ParticleSystem<Vec3>;
template class ParticleSystem<Vec2>;
template class ParticleSystem<Vec1>;

#endif
