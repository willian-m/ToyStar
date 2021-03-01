#include "integrator_rk4.h"
#define BOOST_TEST_MODULE RK4Test
#include <boost/test/unit_test.hpp>

#include <vector>
#include <stdio.h>
#include <iomanip> 
#include <fstream>

#include "vec3.h"
#include "eos_polytropic.h"
#include "particle_system.h"



double get_kinetic_energy(ParticleSystem<Vec3>* sys){
    int npart =  sys->get_nparticles();
    double E = 0;
    for (int ipart=0; ipart<npart;++ipart){
        Particle<Vec3>* particle = sys->get_particle(ipart);
        Vec3 v = particle->get_velocity();
        E += v.x*v.x + v.y*v.y + v.z*v.z;
    }
    return E;
}

double get_potential_energy(ParticleSystem<Vec3>* sys){
    int npart =  sys->get_nparticles();
    double E = 0;
    for (int ipart=0; ipart<npart;++ipart){
        Particle<Vec3>* particle = sys->get_particle(ipart);
        Vec3 r = particle->get_position();
        E += r.x*r.x + r.y*r.y + r.z*r.z;
    }
    return E;
}

double get_total_energy(ParticleSystem<Vec3>* sys){
    int npart =  sys->get_nparticles();
    double E = 0;
    for (int ipart=0; ipart<npart;++ipart){
        Particle<Vec3>* particle = sys->get_particle(ipart);
        Vec3 v = particle->get_velocity();
        Vec3 r = particle->get_position();
        E += v.x*v.x + v.y*v.y + v.z*v.z + r.x*r.x + r.y*r.y + r.z*r.z;
    }
    return E;
}


BOOST_AUTO_TEST_CASE( harmonic_oscillator_test ){

    //setup particle system
    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3> r;
    std::vector<Vec3> v;
    std::vector<double> m;

    Vec3 pos(1., .0, .0);
    Vec3 vel(.0, .0, .0);

    r.push_back(pos);
    v.push_back(vel);
    m.push_back(1.0);


    ParticleSystem<Vec3> current(r,v,m,.1,1.,.0,eos);
    ParticleSystem<Vec3> next(r,v,m,.1,1.,.0,eos);
    ParticleSystem<Vec3> buffer(r,v,m,.1,1.,.0,eos);

    //Integrator parameters
    double total_time = 2.0*M_PI;
    int nsteps = 100000;
    double dt = total_time/nsteps;
    
    //Init integrator
    IntegratorRK4<Vec3> integrator(&current,&next,&buffer,dt);

    std::ofstream total_E_file("Energy_HO_test.txt",std::ios::trunc);
    std::ofstream total_K_file("Kinetic_HO_test.txt",std::ios::trunc);
    std::ofstream total_pot_file("Potential_HO_test.txt",std::ios::trunc);
    total_E_file << "#t\tE" << std::endl;
    total_K_file << "#t\tE" << std::endl;
    total_pot_file << "#t\tE" << std::endl;

    double E0 = get_total_energy(&current);
    total_E_file << "0\t" << E0 << std::endl;
    total_K_file << "0\t" << get_kinetic_energy(&current) << std::endl;
    total_pot_file << "0\t" << get_potential_energy(&current) << std::endl;


    


    for(int istep=0; istep<nsteps;++istep){
        double E1 = get_total_energy(&current);
        integrator.do_step();
        BOOST_CHECK_MESSAGE( abs(E1 - E0) < TOL*100, "Energy changed more than tolerable");
        total_E_file << istep*dt + dt <<"\t"<< E1 << std::endl;
        total_K_file << istep*dt + dt <<"\t" << get_kinetic_energy(&current) << std::endl;
        total_pot_file << istep*dt + dt <<"\t" << get_potential_energy(&current) << std::endl;
        if ((istep == nsteps/4) || (istep == 3*nsteps/4)){
            BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_position().x ) < TOL*100,
                        "Particle at position: "<<current.get_particle(0)-> get_position().x<<
                        ". Expected at 0.");
        } else if (istep == nsteps/2){
            BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_position().x + 1 ) < TOL*100,
                        "Particle at position: "<<current.get_particle(0)-> get_position().x<<
                        ". Expected at -1.");
            BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_velocity().x) < TOL*100,
                        "Particle with velocity: "<<current.get_particle(0)-> get_position().x<<
                        ". Expected at 0.");
        }
    }
    total_E_file << std::flush;
    total_E_file.close();
    total_pot_file << std::flush;
    total_pot_file.close();
    total_K_file << std::flush;
    total_K_file.close();
    std::cout << "Position deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)-> get_position().x - 1.) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_position().x - 1.) < TOL*100,
                        "Particle at position: "<<current.get_particle(0)-> get_position().x<<
                        ". Expected at 1.");

    std::cout << "Velocity deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)-> get_velocity().x) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_velocity().x) < TOL*100,
                        "Particle with velocity: "<<current.get_particle(0)-> get_velocity().x<<
                        ". Expected it to be at rest.");

    
}

BOOST_AUTO_TEST_CASE( damped_harmonic_oscillator ){

    //setup particle system
    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3> r;
    std::vector<Vec3> v;
    std::vector<double> m;

    Vec3 pos(2., .0, .0);
    Vec3 vel(3., .0, .0);

    r.push_back(pos);
    v.push_back(vel);
    m.push_back(1.0);


    ParticleSystem<Vec3> current(r,v,m,.1,1.,1.0,eos);
    ParticleSystem<Vec3> next(r,v,m,.1,1.,1.0,eos);
    ParticleSystem<Vec3> buffer(r,v,m,.1,1.,1.0,eos);

    //Integrator parameters
    double total_time = 2.0*M_PI*20.;
    int nsteps = 10000;
    double dt = total_time/nsteps;
    
    //Init integrator
    IntegratorRK4<Vec3> integrator(&current,&next,&buffer,dt);
    
    
    for(int istep=0; istep<nsteps;++istep){
        integrator.do_step();
    }
    std::cout << std::endl;
    std::cout << "Position deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)-> get_position().x ) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_position().x ) < TOL,
                        "Particle at position: "<<current.get_particle(0)-> get_position().x<<
                        ". Expected at 0.");

        std::cout << "Velocity deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)-> get_velocity().x) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_velocity().x) < TOL,
                        "Particle with velocity: "<<current.get_particle(0)-> get_velocity().x<<
                        ". Expected it to be at rest.");

    
}

/*BOOST_AUTO_TEST_CASE( update_system_test ){
    //setup particle system
    EOSPolytropic* eos = new EOSPolytropic(0.,1.);

    int const n_points = 10; //Test will generate this many random points.
                             //The distance between them should always be positive
    std::vector<Vec3> pos1, vel1;
    std::vector<Vec3> pos2, vel2;
    std::vector<double> mass1, mass2;
    Vec3 pos, vel;


    std::srand(8270509);

    for (int ipoint=0; ipoint<n_points;++ipoint){
        pos.x = (double)(rand()%100);
        pos.y = (double)(rand()%100);
        pos.z = (double)(rand()%100);
        pos1.push_back(pos);
        vel.x = (double)(rand()%100);
        vel.y = (double)(rand()%100);
        vel.z = (double)(rand()%100);
        vel1.push_back(vel);
        double mass = (double)(rand()%100);
        mass1.push_back(mass);
        pos.x = (double)(rand()%100);
        pos.y = (double)(rand()%100);
        pos.z = (double)(rand()%100);
        pos2.push_back(pos);
        vel.x = (double)(rand()%100);
        vel.y = (double)(rand()%100);
        vel.z = (double)(rand()%100);
        vel2.push_back(vel);
        mass = (double)(rand()%100);
        mass2.push_back(mass);
    }

    ParticleSystem<Vec3> current(pos1,vel1,mass1,.05,1.,0.,eos);
    ParticleSystem<Vec3> next(pos2,vel2,mass2,.05,1.,0.,eos);
    ParticleSystem<Vec3> buffer(pos2,vel2,mass2,.05,1.,0.,eos);

    std::cout << &current <<"," <<&next<<std::endl;
    

    IntegratorRK4 integrator(&current,&next,&buffer,1.);
    integrator.update_system();
    std::cout << &current <<"," <<&next<<std::endl;
    
    //Check if we swapped the systems
    for (int ipart=0;ipart<n_points;++ipart){
        Particle<Vec3>* next_part = next.get_particle(ipart);
        BOOST_CHECK_MESSAGE( next_part-> get_position().x == pos1[ipart].x,
            "Wrong particle position in x direction. Expected "<<pos1[ipart].x<<". Received: "<<next_part-> get_position().x);
        BOOST_CHECK_MESSAGE( next_part->get_y() == pos1[ipart].y,
            "Wrong particle position in y direction. Expected "<<pos1[ipart].y<<". Received: "<<next_part->get_y());
        BOOST_CHECK_MESSAGE( next_part->get_z() == pos1[ipart].z,
            "Wrong particle position in z direction. Expected "<<pos1[ipart].z<<". Received: "<<next_part->get_z());
    }
    
}*/