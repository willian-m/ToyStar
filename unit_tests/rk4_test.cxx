#include "integrator_rk4.h"
#define BOOST_TEST_MODULE RK4Test
#include <boost/test/unit_test.hpp>

#include <vector>
#include <stdio.h>
#include <iomanip> 

#include "vec3.h"
#include "eos_polytropic.h"
#include "particle_system.h"

BOOST_AUTO_TEST_CASE( harmonic_oscillator_test ){

    //setup particle system
    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;

    Vec3<double> pos(1., .0, .0);
    Vec3<double> vel(.0, .0, .0);

    r.push_back(pos);
    v.push_back(vel);
    m.push_back(1.0);


    ParticleSystem current(r,v,m,.1,1.,.0,eos);
    ParticleSystem next(r,v,m,.1,1.,.0,eos);
    ParticleSystem buffer(r,v,m,.1,1.,.0,eos);

    //Integrator parameters
    double total_time = 2.0*M_PI;
    int nsteps = 100000;
    double dt = total_time/nsteps;
    
    //Init integrator
    IntegratorRK4 integrator(&current,&next,&buffer,dt);
    

    for(int istep=0; istep<nsteps;++istep){
        integrator.do_step();
        integrator.update_system();
    }

    std::cout << "Position deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)->get_x() - 1.) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_x() - 1.) < TOL*100,
                        "Particle at position: "<<current.get_particle(0)->get_x()<<
                        ". Expected at 1.");

    std::cout << "Velocity deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)->get_vx()) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_vx()) < TOL*100,
                        "Particle with velocity: "<<current.get_particle(0)->get_vx()<<
                        ". Expected it to be at rest.");

    
}

BOOST_AUTO_TEST_CASE( damped_harmonic_oscillator ){

    //setup particle system
    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;

    Vec3<double> pos(2., .0, .0);
    Vec3<double> vel(3., .0, .0);

    r.push_back(pos);
    v.push_back(vel);
    m.push_back(1.0);


    ParticleSystem current(r,v,m,.1,1.,1.0,eos);
    ParticleSystem next(r,v,m,.1,1.,1.0,eos);
    ParticleSystem buffer(r,v,m,.1,1.,1.0,eos);

    //Integrator parameters
    double total_time = 2.0*M_PI*20.;
    int nsteps = 10000;
    double dt = total_time/nsteps;
    
    //Init integrator
    IntegratorRK4 integrator(&current,&next,&buffer,dt);
    
    
    for(int istep=0; istep<nsteps;++istep){
        integrator.do_step();
        integrator.update_system();
    }
    std::cout << std::endl;
    std::cout << "Position deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)->get_x() ) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_x() ) < TOL,
                        "Particle at position: "<<current.get_particle(0)->get_x()<<
                        ". Expected at 0.");

        std::cout << "Velocity deviation from expected: "<< std::setprecision(15) << abs(current.get_particle(0)->get_vx()) << std::endl;
    BOOST_CHECK_MESSAGE( abs(current.get_particle(0)->get_vx()) < TOL,
                        "Particle with velocity: "<<current.get_particle(0)->get_vx()<<
                        ". Expected it to be at rest.");

    
}