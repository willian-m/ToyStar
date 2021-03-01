#include "sph_math.h"
#define BOOST_TEST_MODULE SPHMathTest
#include <boost/test/unit_test.hpp>

#include<stdlib.h>
#include<array>
#include <iostream>
#include <iomanip> 

#include "vec3.h"

#define TOL 1.E-15

BOOST_AUTO_TEST_CASE( positive_distance_test ){

    int const n_points = 10; //Test will generate this many random points.
                             //The distance between them should always be positive
    std::array<Vec3,n_points> pos;

    std::srand(8270509);

    for (int ipoint=0; ipoint<n_points;++ipoint){
        pos[ipoint].x = (double)(rand()%100);
        pos[ipoint].y = (double)(rand()%100);
        pos[ipoint].z = (double)(rand()%100);
    }

    for (int ipoint=0; ipoint<n_points;++ipoint)
        for (int jpoint=ipoint+1; jpoint<n_points;++jpoint)
            BOOST_CHECK( SPHMath::distance(pos[ipoint], pos[jpoint]) > 0. );
        

};

BOOST_AUTO_TEST_CASE( distance_test ){

   Vec3 pos1(3.,.0,.0);
   Vec3 pos2(0,4.,.0);

   BOOST_CHECK_MESSAGE( std::abs(SPHMath::distance(pos1, pos2)-5.) < TOL,
                                 "Distance is: " << SPHMath::distance(pos1, pos2)<<
                                 ". Expected value is 5.");

   pos2.y = 0.;
   pos2.z = 4.;

   BOOST_CHECK_MESSAGE( std::abs(SPHMath::distance(pos1, pos2)-5.) < TOL,
                                 "Distance is: " << SPHMath::distance(pos1, pos2)<<
                                 ". Expected value is 5.");

}

BOOST_AUTO_TEST_CASE( kernel_spline3D_test ){

    Vec3 ri(.3,.0,.0);
    Vec3 rj(.0,.4,.0);
    double h =1.;

    BOOST_CHECK_MESSAGE ( abs(SPHMath::kernel_spline(ri, rj, h) - 2.287852306945996e-01) < TOL,
                          "Spline Value is: " << SPHMath::kernel_spline(ri, rj, h)<<". Expected 2.287852306945996e-01");
    
    ri.x=1.;
    rj.y=1.;

    BOOST_CHECK_MESSAGE ( abs(SPHMath::kernel_spline(ri, rj, h) - 1.599587764401773e-02) < TOL,
                          "Spline Value is: " << SPHMath::kernel_spline(ri, rj, h)<<". Expected 1.599587764401773e-02");

    ri.x=4.;
    BOOST_CHECK_MESSAGE ( abs(SPHMath::kernel_spline(ri, rj, h) ) < TOL,
                          "Spline Value is: " << SPHMath::kernel_spline(ri, rj, h)<<". Expected 0.0");
    
}

BOOST_AUTO_TEST_CASE( kernel_spline3D_norm_test ){

    Vec3 ri(.3,.0,.0);
    Vec3 r_zero(.0,.0,.0);

    double step=8.e-3;
    double L=4;
    double n = L/step-1;

    double sum = .0;
    double h = 1.;
    
    for(int ix=0;ix<n;++ix){
        ri.x = (ix+.5)*step - L/2.;
        for(int iy=0;iy<n;++iy){
            ri.y = (iy+.5)*step - L/2.;
            for(int iz=0;iz<n;++iz){
                ri.z = (iz+.5)*step - L/2.;
                sum += SPHMath::kernel_spline(ri, r_zero, h);
            }
        }
    }
    BOOST_CHECK_MESSAGE( std::abs(sum*step*step*step - 1.) < TOL*1.E+4, //There is an integration error as well, thus we are more lax here
                         "Interpolator does not normalize to 1. Value: "
                         << std::fixed << std::setprecision(15)<< sum*step*step*step);
}

BOOST_AUTO_TEST_CASE( gradient_kernel_spline3D ){

    Vec3 ri(.0,.0,.0);
    Vec3 rj(1.,.0,.0);

    Vec3 spline_grad = SPHMath::gradient_kernel_spline(ri,rj,1);

    BOOST_CHECK_MESSAGE( abs(spline_grad.x  + 3./(4.*M_PI)) < TOL,
                    "Wrong value for x component of spline gradient. Expected "<< -3./(4.*M_PI)
                    << ". Got " << spline_grad.x );
    
    rj.x = .0; rj.y = 1.;
    spline_grad = SPHMath::gradient_kernel_spline(ri,rj,1);
    BOOST_CHECK_MESSAGE( abs(spline_grad.y  + 3./(4.*M_PI)) < TOL,
                    "Wrong value for y component of spline gradient. Expected "<< -3./(4.*M_PI)
                    << ". Got " << spline_grad.y );

    rj.y = .0; rj.z = 1.;
    spline_grad = SPHMath::gradient_kernel_spline(ri,rj,1);
    BOOST_CHECK_MESSAGE( abs(spline_grad.z  + 3./(4.*M_PI)) < TOL,
                    "Wrong value for z component of spline gradient. Expected "<< -3./(4.*M_PI)
                    << ". Got " << spline_grad.z );
    
}