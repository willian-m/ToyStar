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

template ParticleSystem<Vec3>::ParticleSystem(std::vector<Vec3> r,
                               std::vector<Vec3> v,
                               std::vector<double> mass, double lh,
                               double llambda, double lnu,
                               EOSBase* leos
                               );
template ParticleSystem<Vec2>::ParticleSystem(std::vector<Vec2> r,
                               std::vector<Vec2> v,
                               std::vector<double> mass, double lh,
                               double llambda, double lnu,
                               EOSBase* leos
                               );
template ParticleSystem<Vec1>::ParticleSystem(std::vector<Vec1> r,
                               std::vector<Vec1> v,
                               std::vector<double> mass, double lh,
                               double llambda, double lnu,
                               EOSBase* leos
                               );

template <class T>
void ParticleSystem<T>::update_acceleration(){

    clear_grid();

    for (int ipart=0; ipart<nparticles;++ipart){
        T ri = sph_particles[ipart].get_position();
        T vi = sph_particles[ipart].get_velocity();
        T acc = (-1.)*lambda*ri+(-1.)*nu*vi;
        sph_particles[ipart].set_acceleration(acc);
        sph_particles[ipart].clear_neighbor_list(); //Clear list of neighbors, just in case
    }

    //Setup grid
    //TODO: Implement Interface to decide grid boundaries
    setup_grid(Vec3(-1,-1,-1),Vec3(1,1,1));

    find_neighbor_particles(); //Create l
    update_densities(); //Update the densities of all particles
                        // This does not consume neighbour list


    for (int ipart=0; ipart<nparticles;++ipart){
        
        T ri = sph_particles[ipart].get_position();
        T acc_i = sph_particles[ipart].get_acceleration();

        double rho_i = sph_particles[ipart].get_density();
        double P_i = eos->get_pressure(rho_i);
        double p_over_rho2_i = P_i/(rho_i*rho_i);

        auto begin_neighbor_list = sph_particles[ipart].get_begin_neighbor_list();
        auto end_neighbor_list = sph_particles[ipart].get_end_neighbor_list();
        for (auto neigh_part_it = begin_neighbor_list; 
                  neigh_part_it != end_neighbor_list; ++neigh_part_it){
        //for (int jpart=0; jpart<num_neighbors; ++jpart ){
            Particle<T>* neigh_part = neigh_part_it->particle_ptr;
            T rj = neigh_part->get_position();
            T acc_j = neigh_part->get_acceleration();
            double m = neigh_part->get_mass();
            double rho_j = neigh_part->get_density();
            double P_j = eos->get_pressure(rho_j);
            double p_over_rho2_j = P_j/(rho_j*rho_j);
            T gradient = SPHMath::gradient_kernel_spline(ri,rj,h);
            acc_i = acc_i +(-1.)*m*(p_over_rho2_j+p_over_rho2_i)*gradient;
            acc_j = acc_j + m*(p_over_rho2_j+p_over_rho2_i)*gradient;
            neigh_part->set_acceleration(acc_j);
            sph_particles[ipart].set_acceleration(acc_i);
            //Remove ipart from the neighbour list of jpart to avoid double 
            //counting
            neigh_part->erase_neighbor(&sph_particles[ipart]);
        }
          
    }
}

template <class T>
void ParticleSystem<T>::update_densities(){
    for (auto part = sph_particles.begin(); part != sph_particles.end(); ++part){
        double rho = .0;
        for (auto neigh_part = part->get_begin_neighbor_list();
                  neigh_part != part->get_end_neighbor_list(); ++neigh_part){
            rho += neigh_part->particle_ptr->get_mass() * SPHMath::kernel_spline<T>(neigh_part->distance, h);
        }
        //Includes current particle
        rho += part->get_mass()*SPHMath::kernel_spline<T>(.0, h);
        part->set_density(rho);
    }

}

template <class T>
double ParticleSystem<T>::get_point_density(T ri){
    double density = 0;
    for (int jpart =0; jpart<nparticles;++jpart){
        density += sph_particles[jpart].get_mass()
                   *SPHMath::kernel_spline(ri,sph_particles[jpart].get_position(),h);
    }
    return density;
}
template double ParticleSystem<Vec3>::get_point_density(Vec3 ri);
template double ParticleSystem<Vec2>::get_point_density(Vec2 ri);
template double ParticleSystem<Vec1>::get_point_density(Vec1 ri);


template <class T>
void ParticleSystem<T>::setup_grid(Vec3 grid_min, Vec3 grid_max){
    double grid_step = 2*h;
    nx = (int) round((grid_max.x-grid_min.x)/grid_step);
    lim_x[0] = grid_min.x;
    lim_x[1] = (nx-1)*grid_step;
    

    if ((typeid(Vec3) == typeid(T) ) | (typeid(Vec2) == typeid(T) )) {
        ny = (int) round((grid_max.y-grid_min.y)/grid_step);
        lim_y[0] = grid_min.y;
        lim_y[1] = (ny-1)*grid_step;
    } else {
        ny = 1;
        lim_y[0] = .0;
        lim_y[1] = .0;
    }

    if (typeid(Vec3) == typeid(T) ){
        nz = (int) round((grid_max.z-grid_min.z)/grid_step);
        lim_z[0] = grid_min.z;
        lim_z[1] = (nz-1)*grid_step;
    } else {
        nz = 1;
        lim_z[0] = .0;
        lim_z[1] = .0;
    }
    ncells = nx*ny*nz;


    //Init grid with empty cells
    grid = std::vector<Cell<T>>();
    grid.reserve(ncells);
    //Fill it with empty cells
    for (int icell=0; icell<ncells;++icell) grid.emplace_back(Cell<T>(icell));
    build_cell_neighbour_list();
    add_particles_to_grid();
}

template <class T>
void ParticleSystem<T>::clear_grid(){
    for(auto icell = grid.begin(); icell != grid.end(); ++icell){
        icell->clear_cell();
    }
    grid.clear();
}

/// Fills the particles neighbor list 
template <class T>
void ParticleSystem<T>::find_neighbor_particles(){
    for(auto it_cell = grid.begin(); it_cell != grid.end(); ++it_cell ){ //Loop over cells in the grid
        auto ipart = it_cell->particle_list.begin();
        while (ipart != it_cell->particle_list.end() ){ //Loop over particles in the cell
            T ipart_pos = (*ipart)->get_position();
            // 1) Loop over other particles in the cell
            for (Particle<T>* jpart : it_cell->particle_list ){
                //a) Compute distance between ipart and jpart
                double d = SPHMath::distance(ipart_pos,jpart->get_position());
                //b) If distance < 2h, set particles as neighbors
                if (d < 2*h && d > 1.E-14) {
                    (*ipart)->add_neighbour_particle(jpart,d);
                    jpart->add_neighbour_particle((*ipart),d);
                }
            }
            // 2) Loop over particles in neighbour cells
            for (int neighbour_cell_idx : it_cell->neighbour_cell_list){ 
                //      neighbour_cell !=  cell.neighbour_cell_list.end(); ++neighbour_cell ){
                for (Particle<T>* jpart : grid[neighbour_cell_idx].particle_list ){
                    //a) Compute distance between ipart and jpart
                    double d = SPHMath::distance(ipart_pos,jpart->get_position());
                    //b) If distance < 2h, set particles as neighbors
                    if (d < 2*h ) {
                        (*ipart)->add_neighbour_particle(jpart,d);
                        jpart->add_neighbour_particle((*ipart),d);
                    }
                }
            }
            // 3) Remove particle from the cell, to avoid looking into it again
            ipart = it_cell->particle_list.erase(ipart); // This already update to look at next particle
        }
        // TODO: 4) Look at neighbouring cells and remove the current one from its neighbour lists
        // Not critical, since the cell should be empty by now
    }
}

/// Returns the cell index of a point of the grid
template<class T>
int ParticleSystem<T>::get_grid_idx(Vec1 r){
    int ix = (int) ((r.x-lim_x[0])*nx/(lim_x[1]-lim_x[0]));
    ix = ix < 0 ? 0 : (ix >= nx ? nx-1 : ix);

    return ix;
}


template<class T>
int ParticleSystem<T>::get_grid_idx(Vec2 r){
    int ix = (int) ((r.x-lim_x[0])*nx/(lim_x[1]-lim_x[0]));
    ix = ix < 0 ? 0 : (ix >= nx ? nx-1 : ix);

    int iy = (int) ((r.y-lim_y[0])*ny/(lim_y[1]-lim_y[0]));
    iy = iy < 0 ? 0 : (iy >= ny ? ny-1 : iy);
    
    return ix + iy*nx;
}

/// Returns the cell index of a point of the grid
template<class T>
int ParticleSystem<T>::get_grid_idx(Vec3 r){

    int iz = (int) ((r.z-lim_z[0])*nz/(lim_z[1]-lim_z[0]));
    iz = iz < 0 ? 0 : (iz >= nz ? nz-1 : iz);

    return iz*nx*ny + get_grid_idx(Vec2(r.x,r.y));
}

/// Organize the particles inside the grid
template <class T>
void ParticleSystem<T>::add_particles_to_grid(){
    for (auto part = sph_particles.begin(); part != sph_particles.end(); ++part){
        int particle_grid_pos = get_grid_idx( (*part).get_position());
        grid[particle_grid_pos].particle_list.emplace_back(&(*part));
    }
}

/// Fills the neighbour list of the cells
template <class T>
void ParticleSystem<T>::build_cell_neighbour_list(){
    for (int iz = 0; iz<nz; ++iz){
        int iz_plus  = (iz+1 - ((1 + iz)/nz)*nz)*nx*ny;
        int iz_minus = (iz-1 + ((nz - iz)/nz)*nz)*nx*ny;
        int idx_z = iz*nx*ny;
        for (int iy = 0; iy<ny; ++iy){
            int iy_plus  = (iy+1 - ((1 + iy)/ny)*ny)*nx;
            int iy_minus = (iy-1 + ((ny - iy)/ny)*ny)*nx;
            int idx_y = iy*nx;
            for (int ix = 0; ix<nx; ++ix){
                int ix_plus  = ix+1 - ((1 + ix)/nx)*nx;
                int ix_minus = ix-1 + ((nx - ix)/nx)*nx;
                int idx_x = ix;
                int current_cell_idx = idx_z + idx_y + idx_x;
                
                grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + iy_minus + ix_minus );
                if (typeid(T) == typeid(Vec2) )
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + iy_minus + idx_x );
                grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + iy_minus + ix_plus );
                if (typeid(T) == typeid(Vec2) ){
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + idx_y + ix_minus );
                    if (typeid(T) == typeid(Vec3) ) //Avoids self inclusion
                        grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + idx_y + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + idx_y + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + iy_plus + ix_minus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + iy_plus + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_minus + iy_plus + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + iy_minus + ix_minus );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + iy_minus + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + iy_minus + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + idx_y + ix_minus );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + idx_y + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + iy_plus + ix_minus );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + iy_plus + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( idx_z + iy_plus + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + iy_minus + ix_minus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + iy_minus + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + iy_minus + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + idx_y + ix_minus );
                    if (typeid(T) == typeid(Vec3) ) //Avoids self inclusion
                        grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + idx_y + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + idx_y + ix_plus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + iy_plus + ix_minus );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + iy_plus + idx_x );
                    grid[current_cell_idx].neighbour_cell_list.insert( iz_plus + iy_plus + ix_plus );
                }
            }
        }
    }
}
