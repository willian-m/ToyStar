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
template Particle<Vec2>::Particle(Vec2 r_in, Vec2 v_in, Vec2 a_in, double mass);
template Particle<Vec3>::Particle(Vec3 r_in, Vec3 v_in, Vec3 a_in, double mass);
template Particle<Vec1>::Particle(Vec1 r_in, Vec1 v_in, Vec1 a_in, double mass);

template <class T>
void Particle<T>::add_neighbour_particle(Particle<T>* part, double d){
    //NeighbourParticle part_struct = {part, d};
    neighbour_list.emplace_back(NeighbourParticle{part,d});
    
};
template void Particle<Vec1>::add_neighbour_particle(Particle<Vec1>*  part, double d);
template void Particle<Vec2>::add_neighbour_particle(Particle<Vec2>*  part, double d);
template void Particle<Vec3>::add_neighbour_particle(Particle<Vec3>*  part, double d);

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
template void Particle<Vec2>::erase_neighbor(Particle<Vec2>* part);
template void Particle<Vec3>::erase_neighbor(Particle<Vec3>* part);
template void Particle<Vec1>::erase_neighbor(Particle<Vec1>* part);
