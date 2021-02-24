#include "particle_system.h"
#define BOOST_TEST_MODULE ParticleSystemTest
#include <boost/test/unit_test.hpp>

#include <stdio.h>
#include <vector>

#include "vec3.h"
#include "eos_polytropic.h"

#define TOL 1.E-15

BOOST_AUTO_TEST_CASE( density_test ){

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;
    
    //particle 1
    r.push_back(Vec3<double>(.5,0,0));
    v.push_back(Vec3<double>(0,0,0));
    m.push_back(1.);

    //particle 2
    r.push_back(Vec3<double>(-.5,0,0));
    v.push_back(Vec3<double>(0,0,0));
    m.push_back(1.);
    
    double h = 1;

    //Creates systems
    ParticleSystem test_sys(r,v,m,h,1.,1., eos);
    
    Vec3<double> r_probe(.0, .0, .0);
    double density_predicted =  SPHMath::kernel_spline3D(Vec3<double>(.5,0.,0.), Vec3<double>(.0,0.,0.), h)
                              + SPHMath::kernel_spline3D(Vec3<double>(-.5,0.,0.), Vec3<double>(.0,0.,0.), h);

    BOOST_CHECK_MESSAGE( abs(test_sys.get_density(r_probe) -  density_predicted) < TOL,
                         "Wrong density computed. Density is "<<test_sys.get_density(r_probe)
                         << ". It should be "<< density_predicted);

    delete eos;

}

BOOST_AUTO_TEST_CASE( acceleration_test ){

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;
    
    //particle 1
    r.push_back(Vec3<double>(.5,0,0));
    v.push_back(Vec3<double>(0,0,0));
    m.push_back(1.);

    //particle 2
    r.push_back(Vec3<double>(-.5,0,0));
    v.push_back(Vec3<double>(0,0,0));
    m.push_back(1.);
    
    double h = 1;

    ParticleSystem test_sys(r,v,m,h,1.,1., eos);
    
    
    double density_predicted = SPHMath::kernel_spline3D(Vec3<double>(0.,0.,0.), Vec3<double>(0.,0.,0.), h);

    test_sys.update_acceleration();

    Vec3<double> acc = test_sys.get_particle(0)->get_acceleration();

    BOOST_CHECK_MESSAGE( abs(acc.y) < TOL, "Non-zero acceleration on y commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.z) < TOL, "Non-zero acceleration on z commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.x + .5 + 6./(4.*M_PI)) < TOL, 
                        "Wrong acceleration on x component. Got "<< acc.x
                        <<". It should be: " << -.5-3./(4.*M_PI) );
    

    delete eos;

}

BOOST_AUTO_TEST_CASE( acceleration_test_2 ){

    EOSPolytropic* eos = new EOSPolytropic(1.,1.);

    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;
    
    //particle 1
    r.push_back(Vec3<double>(.5,0,0));
    v.push_back(Vec3<double>(1.,0,0));
    m.push_back(1.);

    //particle 2
    r.push_back(Vec3<double>(-.5,0,0));
    v.push_back(Vec3<double>(.0,0,0));
    m.push_back(1.);
    
    double h = 1;

    ParticleSystem test_sys(r,v,m,h,1.,1., eos);
    
    
    double density_predicted = SPHMath::kernel_spline3D(Vec3<double>(0.,0.,0.), Vec3<double>(0.,0.,0.), h);

    test_sys.update_acceleration();
    Vec3<double> acc = test_sys.get_particle(0)->get_acceleration();
    

    BOOST_CHECK_MESSAGE( abs(acc.y) < TOL, "Non-zero acceleration on y commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.z) < TOL, "Non-zero acceleration on z commponent" );
    BOOST_CHECK_MESSAGE( abs(acc.x + .5 + 6./(4.*M_PI)+1) < TOL, 
                        "Wrong acceleration on x component. Got "<< acc.x
                        <<". It should be: " << -.5-3./(4.*M_PI) );
    

    delete eos;

}