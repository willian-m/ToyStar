#ifndef CELL_H
#define CELL_H

#include <vector>
#include <array>
#include <set>

#include "vec2.h"
#include "vec3.h"
#include "particle.h"

template <class T>
class Cell{

private:

    int id;
    
public:

    std::vector<Particle<T>*> particle_list;
    std::set<int> neighbour_cell_list;

    Cell(int id);
    void clear_cell();
    ~Cell();


};

template<class T>
inline Cell<T>::Cell(int id):id(id){};
template<class T>
inline Cell<T>::~Cell(){
    particle_list.clear();
    neighbour_cell_list.clear();
};
template <class T>
inline void Cell<T>::clear_cell(){
    particle_list.clear();
    neighbour_cell_list.clear();
}

template class Cell<Vec3>;
template class Cell<Vec2>;
template class Cell<Vec1>;

#endif

