#include "particle_system.h"
#define BOOST_TEST_MODULE ParticleSystemTest
#include <boost/test/unit_test.hpp>

#include <stdio.h>
#include <vector>
#include "vec3.h"

#define TOL 1.E-15

BOOST_AUTO_TEST_CASE( density_test ){

    int nx = 100;
    int ny = 100;

    double Lx = 4;
    double Ly = 4;

    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;
    double h = Lx/nx/10;

    for (int iy = 0; iy < ny; iy++){
        double y = iy*Ly/ny - Ly/2.;
        for (int ix = 0; ix < ny; ix++){
            double x = ix*Lx/nx - Lx/2.;
            r.push_back(Vec3<double>(x,y,0));
            v.push_back(Vec3<double>(0,0,0));
            m.push_back(1.);
        }
    }

    std::cout << "Creatimg system";
    ParticleSystem test_sys(r,v,m,h);
    
    Vec3<double> r_probe(50*Lx/nx - Lx/2., 50*Ly/ny - Ly/2., 0);
    double density_predicted = SPHMath::kernel_spline3D(Vec3<double>(0.,0.,0.), Vec3<double>(0.,0.,0.), h);

    BOOST_CHECK_MESSAGE( abs(test_sys.get_density(r_probe) -  density_predicted) < TOL,
                         "Wrong density computed. Density is "<<test_sys.get_density(r_probe)
                         << ". It should be "<< 1./(4.*M_PI*std::pow(Lx/nx,3)));

}