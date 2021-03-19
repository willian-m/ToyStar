#ifndef CELL_H
#define CELL_H

#include <vector>
#include <array>
#include <set>


#include "particle.h"

template <class T>
class Cell{

private:

    int id;
    
    /*double grid_step;
    int nx, ny, nz;

    std::array<int,2> lim_x, lim_y, lim_z;
    
    std::vector<std::vector<Particle<T>*>> particle_list; //< Array of list of pointers
    std::vector<std::array<int,27>> cell_neighbour_list; //< List of cells in the neighborhood of 

    int get_grid_idx(Vec2 r);
    int get_grid_idx(Vec3 r);
    void add_particle(Particle<T>* particle);
    void build_cell_neighbour_list();*/

public:

    std::vector<Particle<T>*> particle_list;
    std::set<int> neighbour_cell_list;

    Cell(int id);
    void clear_cell();
    ~Cell();
    //void add_neighbour_cell( Cell<T>* neighbour_cell );
    //void add_particle( Particle<T>* particle );
    //const std::vector<Particle<T>*>& get_particle_list();
    //const std::set<Cell<T>*>& get_cell_neighbour_list();

    /*Grid(Vec3 grid_min, Vec3 grid_max, double h);
    Grid(Vec2 grid_min, Vec2 grid_max, double h);
    
    void fill_grid(ParticleSystem<T>* part_sys);
    int get_ncells();
    std::vector<Particle<T>*> get_particle_list(int icell);
    std::vector<Particle<T>*> get_particle_list_neighbour(int icell, int ineigh);*/



};

template<class T>
inline Cell<T>::Cell(int id):id(id){};
template<class T>
inline Cell<T>::~Cell(){
    particle_list.clear();
    neighbour_cell_list.clear();
};
//template<class T>
//inline void Cell<T>::add_neighbour_cell(Cell<T>* neighbour_cell){neighbour_cell_list.insert(neighbour_cell);};
//template<class T>
//inline void Cell<T>::add_particle(Particle<T>* particle){particle_list.emplace_back(particle);};
//template<class T>
//inline const std::vector<Particle<T>*>& Cell<T>::get_particle_list(){return particle_list;};
//template<class T>
//inline const std::set<Cell<T>*>& Cell<T>::get_cell_neighbour_list(){return neighbour_cell_list;};


/*
template<class T>
inline void Grid<T>::fill_grid(ParticleSystem<T>* part_sys){
    for (Particle<T> part : part_sys->sph_particles) add_particle(&part);};

template <class T>
inline int Grid<T>::get_ncells(){return nx*ny*nz;};

template <class T>
inline std::vector<Particle<T>*> Grid<T>::get_particle_list(int icell){ return &particle_list[icell]};*/

template class Cell<Vec3>;
template class Cell<Vec2>;

#endif
