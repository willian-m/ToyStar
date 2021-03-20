#include "particle.h"

///TODO: Add functions documentations

/** \class Particle
 A class to store SPH particle properties. There is no dynamic in here
*/

template <class T>
Particle<T>::Particle(T r_in, T v_in, T a_in, double mass){
    r = r_in;
    v = v_in;
    a = a_in;
    this->mass = mass;
}

template <class T>
void Particle<T>::add_neighbour_particle(Particle<T>* part, double d){
    //NeighbourParticle part_struct = {part, d};
    neighbour_list.emplace_back(NeighbourParticle{part,d});
    
};

template <class T>
void Particle<T>::clear_neighbor_list(){
    neighbour_list.clear();
}

template <class T>
void Particle<T>::erase_neighbor(Particle<T>* part){
    for (auto neigh_part = neighbour_list.begin();
              neigh_part != neighbour_list.end(); ){
        if ( neigh_part->particle_ptr == part ) {
            neighbour_list.erase(neigh_part);
            break;
        } else {
            ++neigh_part;
        }
    }
}