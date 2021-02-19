#include "sph_math.h"

//Norms
const double SPHMath::norm3D = 1./(4.*M_PI);

SPHMath::SPHMath(){};

double SPHMath::kernel_spline3D(Vec3<double> ri, Vec3<double> rj, double h){
    double q = SPHMath::distance(ri, rj)/h;
    double norm = norm3D/std::pow(h,3);


    if (q < 1.){
        return (std::pow(2.-q,3) - 4.*pow(1.-q,3))*norm;
    } else if (q < 2){
        return (std::pow(2.-q,3))*norm;
    } else {
        return 0;
    }

    return 0;

}

double SPHMath::distance(Vec3<double> ri, Vec3<double> rj){
    return(std::sqrt( (ri.x - rj.x)*(ri.x - rj.x)
                     +(ri.y - rj.y)*(ri.y - rj.y)
                     +(ri.z - rj.z)*(ri.z - rj.z) ) );
}
