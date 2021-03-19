#include "particle_system.h"
#define BOOST_TEST_MODULE ParticleSystemTest
#include <boost/test/unit_test.hpp>

#include <stdio.h>
#include <vector>

#include "vec3.h"
#include "eos_polytropic.h"

#define TOL 1.E-15

BOOST_AUTO_TEST_CASE( neighbour_test ){

    //A) Create a 3D system with 400 particles placed randomly
    std::vector<Vec3> r;
    std::vector<Vec3> v;
    std::vector<double> m;

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);


    std::srand(8270509);

    for (int ipoint=0; ipoint<400;++ipoint){
        double x = (double)((rand()-RAND_MAX/2))/RAND_MAX;
        double y = (double)((rand()-RAND_MAX/2))/RAND_MAX;
        double z = (double)((rand()-RAND_MAX/2))/RAND_MAX;

        r.push_back(Vec3(x,y,z));
        v.push_back(Vec3(0,0,0));
        m.push_back(1.);
    }
    double h = 0.05;
    ParticleSystem<Vec3> sys(r, v, m, h, 0.,1.,eos);

    sys.update_acceleration(); //Triggers algorithm of filling particle neighbors

    
    for(int ipart=0; ipart < 400; ++ipart){
        Particle<Vec3>* part = sys.get_particle(ipart);
        std::vector<Particle<Vec3>*> neigh_part;
        //Loop over all particles, finding the ones that are in the neighborhood of ipart
        for(int jpart=0; jpart < 400; ++jpart){
            Particle<Vec3>* part_candidate = sys.get_particle(jpart);
            double d = SPHMath::distance(part->get_position(),part_candidate->get_position());
            if ( d < 2*h && d > 1.e-14 ) neigh_part.push_back(part_candidate);
        }
        //Now we test if all neighbour particles were included in the neighbour list of ipart
        int num_neighbors = part->get_num_neighbors();
        BOOST_CHECK_MESSAGE( num_neighbors == neigh_part.size(),
                            "Particle neighbour list has different size of actual neighbors. "
                          <<"Number of neighbors in list: " << num_neighbors <<". Actual number of neighbors: " <<neigh_part.size());
        for (int ipart_neigh=0; ipart_neigh << num_neighbors; ++ipart_neigh){
            Particle<Vec3>* part_neigh_candidate = part->get_neighbor(ipart_neigh);
            bool particle_found = false;
            for( Particle<Vec3>* actual_part_neigh : neigh_part){
                if (actual_part_neigh == part_neigh_candidate) particle_found = true;
            }
            BOOST_CHECK_MESSAGE(particle_found,"Missing particle in the neighborhood list found.");
        }
         
    }

}

BOOST_AUTO_TEST_CASE( density_test ){

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3> r;
    std::vector<Vec3> v;
    std::vector<double> m;
    
    //particle 1
    r.push_back(Vec3(.5,0,0));
    v.push_back(Vec3(0,0,0));
    m.push_back(1.);

    //particle 2
    r.push_back(Vec3(-.5,0,0));
    v.push_back(Vec3(0,0,0));
    m.push_back(1.);
    
    double h = 1;

    //Creates systems
    ParticleSystem<Vec3> test_sys(r,v,m,h,1.,1., eos);
    
    Vec3 r_probe(.0, .0, .0);
    double density_predicted =  SPHMath::kernel_spline(Vec3(.5,0.,0.), Vec3(.0,0.,0.), h)
                              + SPHMath::kernel_spline(Vec3(-.5,0.,0.), Vec3(.0,0.,0.), h);

    BOOST_CHECK_MESSAGE( abs(test_sys.get_density(r_probe) -  density_predicted) < TOL,
                         "Wrong density computed. Density is "<<test_sys.get_density(r_probe)
                         << ". It should be "<< density_predicted);

    delete eos;

}

BOOST_AUTO_TEST_CASE( acceleration_test ){

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3> r;
    std::vector<Vec3> v;
    std::vector<double> m;
    
    //particle 1
    r.push_back(Vec3(.5,0,0));
    v.push_back(Vec3(0,0,0));
    m.push_back(1.);

    //particle 2
    r.push_back(Vec3(-.5,0,0));
    v.push_back(Vec3(0,0,0));
    m.push_back(1.);
    
    double h = 1;

    ParticleSystem<Vec3> test_sys(r,v,m,h,1.,1., eos);
    
    
    double density_predicted = SPHMath::kernel_spline(Vec3(0.,0.,0.), Vec3(0.,0.,0.), h);

    test_sys.update_acceleration();

    Vec3 acc = test_sys.get_particle(0)->get_acceleration();

    BOOST_CHECK_MESSAGE( abs(acc.y) < TOL, "Non-zero acceleration on y commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.z) < TOL, "Non-zero acceleration on z commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.x + .5 - 6./(4.*M_PI)) < TOL, 
                        "Wrong acceleration on x component. Got "<< acc.x
                        <<". It should be: " << -.5-3./(4.*M_PI) );
    

    delete eos;

}

BOOST_AUTO_TEST_CASE( acceleration_test_2 ){

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3> r;
    std::vector<Vec3> v;
    std::vector<double> m;
    
    //particle 1
    r.push_back(Vec3(.5,0,0));
    v.push_back(Vec3(1.,0,0));
    m.push_back(1.);

    //particle 2
    r.push_back(Vec3(-.5,0,0));
    v.push_back(Vec3(.0,0,0));
    m.push_back(1.);
    
    double h = 1;

    ParticleSystem<Vec3> test_sys(r,v,m,h,1.,1., eos);
    
    
    double density_predicted = SPHMath::kernel_spline(Vec3(0.,0.,0.), Vec3(0.,0.,0.), h);

    test_sys.update_acceleration();
    Vec3 acc = test_sys.get_particle(0)->get_acceleration();
    

    BOOST_CHECK_MESSAGE( abs(acc.y) < TOL, "Non-zero acceleration on y commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.z) < TOL, "Non-zero acceleration on z commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.x + .5 - 6./(4.*M_PI)+1) < TOL, 
                        "Wrong acceleration on x component. Got "<< acc.x
                        <<". It should be: " << -.5-3./(4.*M_PI) );
    

    delete eos;

}