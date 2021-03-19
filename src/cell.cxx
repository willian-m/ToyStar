#include "cell.h"

template <class T>
void Cell<T>::clear_cell(){
    particle_list.clear();
    neighbour_cell_list.clear();
}
/*
template<class T>
Grid<T>::Grid(Vec3 grid_min, Vec3 grid_max, double h):
                grid_step(2.*h){

    nx = (int) round((grid_max.x-grid_min.x)/grid_step);
    ny = (int) round((grid_max.y-grid_min.y)/grid_step);
    nz = (int) round((grid_max.z-grid_min.z)/grid_step);

    lim_x[0] = grid_min.x;
    lim_y[0] = grid_min.y;
    lim_z[0] = grid_min.z;
    
    lim_x[1] = (nx-1)*grid_step;
    lim_y[1] = (ny-1)*grid_step;
    lim_z[1] = (nz-1)*grid_step;

    particle_list = std::vector<std::vector<Particle<T>*>,nx*ny*nz>;
    build_cell_neighbour_list();
}

template<class T>
Grid<T>::Grid(Vec2 grid_min, Vec2 grid_max, double h):
                grid_step(2.*h){
    lim_x[0] = grid_max.x;
    lim_y[0] = grid_max.y;
    
    lim_x[1] = grid_max.x;
    lim_y[1] = grid_max.y;
    
    nx = ((lim_x[1]-lim_x[0])/(grid_step));
    ny = ((lim_y[1]-lim_y[0])/(grid_step));
    nz = 1;

    particle_list = std::vector<std::vector<Particle<T>*>,nx*ny>
    build_cell_neighbour_list();
}

template<class T>
int Grid<T>::get_grid_idx(Vec2 r){
    int ix = (int) (r.x-lim_x[0])/nx;
    ix = ix < 0 ? 0 : (ix >= nz ? nx-1 : ix);

    int iy = (int) (r.y-lim_y[0])/ny;
    iy = iy < 0 ? 0 : (iy >= ny ? ny-1 : iy);
    
    return ix + iy*nx;
}

template<class T>
int Grid<T>::get_grid_idx(Vec3 r){

    int iz = (int) (r.z-lim_z[0])/nz;
    iz = iz < 0 ? 0 : (iz >= nz ? nz-1 : iz);

    return iz*nx*ny + get_grid_idx(Vec2(r.x,r.y));
}

template<class T>
void Grid<T>::add_particle(Particle<T>* particle){
    int particle_grid_pos = get_grid_idx(particle->get_position())
    particle_list[particle_grid_pos].emplace_back(particle);
}

template<class T>
void Grid<T>::build_cell_neighbour_list(){
    for (int iz=0;iz<nz;++iz){
        iz_plus = (iz+1 + ((1 + iz)/nz)*nz)*nx*ny;
        idx_z = iz*nx*ny;
        for (int iy=0;iy<ny;++iy){
            iy_plus = (iy + 1 + ((1 + iy)/ny)*ny)*nx;
            idx_y = iy*nx + iz*nx*ny;
            for (int ix=0;ix<nx;++ix){
                std::array<int,26>
                ix_plus = ix+1 + ((1 + ix)/nx)*nx;
                idx = ix + idx_y;
                //Store neighbors index
                cell_neighbour_list.emplace_back(std::array<int,26>);
                cell_neighbour_list[idx][0] = idx_z + idx_y + ix;
                cell_neighbour_list[idx][1] = idx_z + idx_y + ix_plus;
                cell_neighbour_list[idx][2] = idx_z + iy_plus + ix;
                cell_neighbour_list[idx][3] = idx_z + iy_plus + ix_plus;
                cell_neighbour_list[idx][4] = iz_plus + idx_y + ix;
                cell_neighbour_list[idx][5] = iz_plus + idx_y + ix_plus;
                cell_neighbour_list[idx][6] = iz_plus + iy_plus + ix;
                cell_neighbour_list[idx][7] = iz_plus + iy_plus + ix_plus;   
            }                
        }
    }
}

template<class T>
std::vector<Particle<T>*> Grid<T>::get_particle_list_neighbour(int icell, int ineigh){
    return get_particle_list(cell_neighbour_list[icell][ineigh]);
}
*/