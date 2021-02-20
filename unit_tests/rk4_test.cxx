#include "integrator_rk4.h"
#define BOOST_TEST_MODULE RK4Test
#include <boost/test/unit_test.hpp>

#include <vector>
#include <stdio.h>

#include "vec3.h"
#include "eos_polytropic.h"
#include "particle_system.h"

#define TOL 1.E-6

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

    //Integrator parameters
    double total_time = 2.0*M_PI;
    int nsteps = 10000000;
    double dt = total_time/nsteps;
    
    //Init integrator
    IntegratorRK4 integrator(&current,&next,dt);
    

    for(int istep; istep<nsteps;++istep){
        integrator.do_step();
        integrator.update_system();
    }

    BOOST_CHECK_MESSAGE( abs(next.get_particle(0)->get_x() - 1.) < TOL,
                        "Particle at position: "<<next.get_particle(0)->get_x()<<
                        ". Expected at 1.");

    BOOST_CHECK_MESSAGE( abs(next.get_particle(0)->get_vx()) < TOL,
                        "Particle with velocity: "<<next.get_particle(0)->get_vx()<<
                        ". Expected it to be at rest.");

    
}

