#include "sph_math.h"

//Norms
const double SPHMath::norm3D = 1./(4.*M_PI);
const double SPHMath::norm2D = 5./(14.*M_PI);
SPHMath::SPHMath(){};

double SPHMath::kernel_spline(Vec3 ri, Vec3 rj, double h){
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

template <class T>
double SPHMath::kernel_spline(double d, double h){
    double q = d/h;
    
    double norm = ((typeid(T) == typeid(Vec3) ) ? norm3D/std::pow(h,3) : norm2D/std::pow(h,2));
    
    if (q < 1.){
        return (std::pow(2.-q,3) - 4.*pow(1.-q,3))*norm;
    } else if (q < 2){
        return (std::pow(2.-q,3))*norm;
    } else {
        return 0;
    }

    return 0;
}
template double SPHMath::kernel_spline<Vec3>(double d, double h);
template double SPHMath::kernel_spline<Vec2>(double d, double h);


double SPHMath::kernel_spline(Vec2 ri, Vec2 rj, double h){
    double q = SPHMath::distance(ri, rj)/h;
    double norm = norm2D/std::pow(h,2);


    if (q < 1.){
        return (std::pow(2.-q,3) - 4.*pow(1.-q,3))*norm;
    } else if (q < 2){
        return (std::pow(2.-q,3))*norm;
    } else {
        return 0;
    }

    return 0;
}

double SPHMath::kernel_spline3D(double distance, double h){
    double q = distance/h;
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

double SPHMath::kernel_spline2D(double distance, double h){
    double q = distance/h;
    double norm = norm2D/std::pow(h,2);

    if (q < 1.){
        return (std::pow(2.-q,3) - 4.*pow(1.-q,3))*norm;
    } else if (q < 2){
        return (std::pow(2.-q,3))*norm;
    } else {
        return 0;
    }

    return 0;
}


double SPHMath::distance(Vec3 ri, Vec3 rj){
    return(std::sqrt( (ri.x - rj.x)*(ri.x - rj.x)
                     +(ri.y - rj.y)*(ri.y - rj.y)
                     +(ri.z - rj.z)*(ri.z - rj.z) ) );
}

double SPHMath::distance(Vec2 ri, Vec2 rj){
    return(std::sqrt( (ri.x - rj.x)*(ri.x - rj.x)
                     +(ri.y - rj.y)*(ri.y - rj.y) ) );
}


Vec3 SPHMath::gradient_kernel_spline(Vec3 ri, Vec3 rj, double h){
    double d = SPHMath::distance(ri, rj);

    if (d < TOL) return Vec3(.0,.0,.0); //Avoid division by 0 at the origin

    double q = d/h;
    double norm = norm3D/std::pow(h,4);
    
    double del_W_del_r = .0;


    if (q < 1.){
        del_W_del_r = -(3.*std::pow(2.-q,2) - 12.*pow(1.-q,2))*norm;
    } else if (q < 2.){
        del_W_del_r = -(3.*std::pow(2.-q,2))*norm;
    }
    return Vec3( del_W_del_r*(ri.x-rj.x)/d, del_W_del_r*(ri.y-rj.y)/d,
                         del_W_del_r*(ri.z-rj.z)/d ) ;
}

Vec2 SPHMath::gradient_kernel_spline(Vec2 ri, Vec2 rj, double h){
    double d = SPHMath::distance(ri, rj);

    if (d < TOL) return Vec3(.0,.0,.0); //Avoid division by 0 at the origin

    double q = d/h;
    double norm = norm2D/std::pow(h,3);
    
    double del_W_del_r = .0;


    if (q < 1.){
        del_W_del_r = -(3.*std::pow(2.-q,2) - 12.*pow(1.-q,2))*norm;
    } else if (q < 2.){
        del_W_del_r = -(3.*std::pow(2.-q,2))*norm;
    }
    return Vec2( del_W_del_r*(ri.x-rj.x)/d, del_W_del_r*(ri.y-rj.y)/d ) ;
}