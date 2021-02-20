#include <vector>
#include <math.h> //Gamma function
#include <iostream>
#include <fstream>

#include "eos_polytropic.h"
#include "vec3.h"
#include "particle_system.h"
#include "integrator_rk4.h"

int main(int argc, char* argv[]){

    //System parameters
    const double total_time = 16*10;
    const int nsteps = 400*10;
    const double dt = total_time/nsteps;

    //System parameters - initial positions
    const double L = 2.;
    const int nside = 8;
    const double dL = L/nside;

    //System parameters
    const double star_mass = 2.;
    const double star_radius = .75;
    const double particle_mass = star_mass/pow(nside,3);
    const double smoothing = 0.04/sqrt(pow(nside,3)/1000.);
    const double damping = 1.;

    //Allocates eos
    const double pressure_const = 0.1;
    const double poly_const = 1;
    EOSPolytropic* eos = new EOSPolytropic(pressure_const,poly_const);

    //Get attractive constant
    const double fact1 = 2*pressure_const*(1+poly_const)*pow(M_PI,-3./(2*poly_const));
    const double fact2 = tgamma(5./2. + poly_const)*star_mass/(tgamma(1.+poly_const)*pow(star_radius,3));
    double lambda = fact1*fact2/pow(star_mass,2);


    //Creates particles
    std::vector<Vec3<double>> r;
    std::vector<Vec3<double>> v;
    std::vector<double> m;
    Vec3<double> pos(.0, .0, .0);
    Vec3<double> vel(.0, .0, .0);
    for (int ix = 0; ix<nside; ++ix){
        pos.x = -L/2 + ix*dL;
        for (int iy = 0; iy<nside; ++iy){
            pos.y = -L/2 + iy*dL;
            for (int iz = 0; iz<nside; ++iz){
                pos.z = -L/2 + iz*dL;
                r.push_back(pos);
                v.push_back(vel);
                m.push_back(particle_mass);
            }
        }
    }




    //Init particle system
    ParticleSystem current(r,v,m,smoothing,lambda,damping,eos);
    ParticleSystem next(r,v,m,smoothing,lambda,damping,eos);

    //Init integrator
    IntegratorRK4 integrator(&current,&next,dt);
    
    //Simulate it!    
    int barWidth = 70; //For progress bar
    for(int istep=0; istep<nsteps;++istep){
        
        //Evolve the system
        integrator.do_step();
        integrator.update_system();

        if (istep%10 != 0) continue;
        //Progress bar: https://stackoverflow.com/a/14539953/2754579
        float progress = (double) istep/(double) nsteps;
        std::cout << "[";
        int bpos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < bpos) std::cout << "=";
            else if (i == bpos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
        
    }
    std::cout << std::endl;
    //Iterate over particles, getting their position, velocity, mass and density
    //Dump it in a text file
    std::ofstream out_file("output.txt",std::ios::trunc);

    if (out_file.is_open() && out_file.good()){
        out_file << "# mass\tx\ty\tz\tvx\tvy\tvz\tdensity"<<std::endl;
        int nparticles = current.get_nparticles();
        for (int ipart=0; ipart<nparticles;++ipart){
            Particle* p = current.get_particle(ipart);
            double x = p->get_x(); double y = p->get_y(); double z = p->get_z();
            out_file <<p->get_mass()<<"\t"
                     <<x<<"\t"<<y<<"\t"<<z<<"\t"
                     <<p->get_vx()<<"\t"<<p->get_vy()<<"\t"<<p->get_vz()<<"\t"
                     <<current.get_density(Vec3<double>(x,y,z))<<std::endl;
        }
    }
    out_file << std::flush;
    out_file.close();
    

    //In addition, gets the density at (z=0, y=0)
    std::ofstream out_density("density.txt",std::ios::trunc);

    if (out_density.is_open() && out_density.good()){
        out_density << "# r\tdensity"<<std::endl;
        double const max_r=2.;
        double const dr=1.e-2;
        int const nr = max_r/dr;
        for (int ir=0; ir<nr;++ir){
            double r = ir*dr;
            out_density <<r<<"\t"<<current.get_density(Vec3<double>(r,0,0))<<std::endl;
        }
    }
    out_density << std::flush;
    out_density.close();

    




}