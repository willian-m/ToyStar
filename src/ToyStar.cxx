#include <vector>
#include <math.h> //Gamma function
#include <iostream>
#include <fstream>

#include "eos_polytropic.h"
#include "vec3.h"
#include "particle_system.h"
#include "integrator_rk4.h"

enum HydroType {
    Bidimensional = 2,
    Tridimensional = 3
};

struct ToyStarPars
{
    double L;
    int nside;
    double dL;

    double star_mass;
    double star_radius;
    double particle_mass;

    double smoothing;
    double damping;

    double pressure_const;
    double poly_const;

    double lambda;

    double total_time;
    double nsteps;
    double dt;

    HydroType mode;
};

void PrintParticleSystem(ParticleSystem<Vec3>* sys,std::string filename){
    std::ofstream out_file(filename,std::ios::trunc);

    //Dumps the data
    if (out_file.is_open() && out_file.good()){
        out_file << "# mass\tx\ty\tz\tvx\tvy\tvz\tdensity"<<std::endl;
        int nparticles = sys->get_nparticles();
        for (int ipart=0; ipart<nparticles;++ipart){
            Particle<Vec3>* p = sys->get_particle(ipart);
            Vec3 pos = p->get_position();
            Vec3 vel = p->get_velocity();
            double x = pos.x; double y = pos.y; double z = pos.z;
            out_file <<p->get_mass()<<"\t"
                <<x<<"\t"<<y<<"\t"<<z<<"\t"
                <<vel.x<<"\t"<<vel.y<<"\t"<<vel.z<<"\t"
                <<sys->get_density(Vec3(x,y,z))<<std::endl;
        }
    }
    out_file << std::flush;
    out_file.close();
}
void PrintParticleSystem(ParticleSystem<Vec2>* sys,std::string filename){
    std::ofstream out_file(filename,std::ios::trunc);

    //Dumps the data
    if (out_file.is_open() && out_file.good()){
        out_file << "# mass\tx\ty\tz\tvx\tvy\tvz\tdensity"<<std::endl;
        int nparticles = sys->get_nparticles();
        for (int ipart=0; ipart<nparticles;++ipart){
            Particle<Vec2>* p = sys->get_particle(ipart);
            Vec2 pos = p->get_position();
            Vec2 vel = p->get_velocity();
            double x = pos.x; double y = pos.y; double z = .0;
            out_file <<p->get_mass()<<"\t"
                <<x<<"\t"<<y<<"\t"<<z<<"\t"
                <<vel.x<<"\t"<<vel.y<<"\t"<<.0<<"\t"
                <<sys->get_density(Vec2(x,y))<<std::endl;
        }
    }
    out_file << std::flush;
    out_file.close();
}

void PrintDensityAlongX(ParticleSystem<Vec2>* sys,std::string filename){

    //In addition, gets the density at (z=0, y=0)
    std::ofstream out_density(filename,std::ios::trunc);

    if (out_density.is_open() && out_density.good()){
        out_density << "# r\tdensity"<<std::endl;
        double const max_r=2.;
        double const dr=1.e-2;
        int const nr = max_r/dr;
        for (int ir=0; ir<nr;++ir){
            double r = ir*dr;
            out_density <<r<<"\t"<<sys->get_density(Vec2(r,0))<<std::endl;
        }
    }
    out_density << std::flush;
    out_density.close();
}

void PrintDensityAlongX(ParticleSystem<Vec3>* sys,std::string filename){

    //In addition, gets the density at (z=0, y=0)
    std::ofstream out_density(filename,std::ios::trunc);

    if (out_density.is_open() && out_density.good()){
        out_density << "# r\tdensity"<<std::endl;
        double const max_r=2.;
        double const dr=1.e-2;
        int const nr = max_r/dr;
        for (int ir=0; ir<nr;++ir){
            double r = ir*dr;
            out_density <<r<<"\t"<<sys->get_density(Vec3(r,0,0))<<std::endl;
        }
    }
    out_density << std::flush;
    out_density.close();
}


template <class T>
void EvolveIt(ToyStarPars pars, std::vector<T>* r, std::vector<T>* v,   std::vector<double>* m){
    EOSPolytropic* eos = new EOSPolytropic(pars.pressure_const,pars.poly_const);
    //Creates particles
   
    //Each particle behaves as an damped HO. We use our knowledge of the 
    //anallytic solution to estimate how much time we must wait to get 
    //to the stationary solution 

    //Print parameters
    std::cout << "MODE: " << pars.mode << "D"<< std::endl;
    std::cout << "SPH parameters" << std::endl;
    std::cout << "---------------------" << std::endl;
    std::cout << "Number of particles..: " << pow(pars.nside,(int) pars.mode) << std::endl;
    std::cout << "Particle mass........: " << pars.particle_mass << std::endl;
    std::cout << "Smoothing (h)........: " << pars.smoothing << std::endl;

    std::cout << "Star paramneters" << std::endl;
    std::cout << "---------------------" << std::endl;
    std::cout << "Damping..............: " << pars.damping << std::endl;
    std::cout << "Gravity..............: " << pars.lambda << std::endl;
    std::cout << "Star mass............: " << pars.star_mass << std::endl;
    std::cout << "Star radius..........: " << pars.star_radius << std::endl;


    std::cout << "EOS parameters" << std::endl;
    std::cout << "---------------------" << std::endl;
    std::cout << "Pressure constant....: " << pars.pressure_const << std::endl;
    std::cout << "Polytropic constant..: " << pars.poly_const << std::endl;

    std::cout << "Simulation parameters" << std::endl;
    std::cout << "---------------------" << std::endl;
    std::cout << "Simulation time......: " << pars.total_time << std::endl;
    std::cout << "Number of steps......: " << pars.nsteps << std::endl;
    std::cout << "Step size............: " << pars.dt << std::endl;

    //Init particle system
    ParticleSystem<T> current(*r,*v,*m,pars.smoothing,pars.lambda,pars.damping,eos);
    ParticleSystem<T> next(*r,*v,*m,pars.smoothing,pars.lambda,pars.damping,eos);
    ParticleSystem<T> buffer(*r,*v,*m,pars.smoothing,pars.lambda,pars.damping,eos);

    //Init integrator
    IntegratorRK4<T> integrator(&current,&next,&buffer,pars.dt);

    //Simulate it!    
    int barWidth = 70; //Width of progress bar
    for(int istep=0; istep<pars.nsteps;++istep){

        //Evolve the system
        integrator.do_step();
        //Iterate over particles, getting their position, velocity, mass and density
        //Dump it in a text file

        if (istep%10 != 0) continue;
        //Create filename
        std::stringstream filename_stream;
        filename_stream << "output_iter_" << istep << ".txt";
        std::string filename;
        filename_stream >> filename;

        //Print output
        PrintParticleSystem(&current, filename);
        
        //Progress bar: https://stackoverflow.com/a/14539953/2754579
        float progress = (double) istep/(double) pars.nsteps;
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
    PrintDensityAlongX(&current,"density.txt");
}

int main(int argc, char* argv[]){

    //Choose run mode
    ToyStarPars pars;
    pars.mode = Bidimensional;
    
    //System parameters - initial positions
    pars.L = 2;
    pars.nside = 20;
    pars.dL = pars.L/pars.nside;

    //System parameters
    pars.star_mass = 2.;
    pars.star_radius = .75;
    pars.particle_mass = pars.star_mass/pow(pars.nside,(int) pars.mode);
    //const double smoothing = 0.04/sqrt(pow(nside,3)/1000.);
    pars.smoothing = pars.dL;//*1.8;//1.5*0.04/sqrt(pow(nside,3)/1000.);
    pars.damping = 2.;

    //Allocates eos
    pars.pressure_const = .1;
    pars.poly_const = 1;

    //Get attractive constant
    pars.lambda=20;
    if (pars.mode == Bidimensional){
        const double fact1 = 2*pars.pressure_const*pow(M_PI,-1./pars.poly_const);
        const double fact2 = pow(pars.star_mass*(1+pars.poly_const)/pow(pars.star_radius,2),1.+1./pars.poly_const);
        pars.lambda = fact1*fact2/pars.star_mass;
    } else if (pars.mode == Tridimensional){
        const double fact1 = 2*pars.pressure_const*(1+pars.poly_const)*pow(M_PI,-3./(2*pars.poly_const));
        const double fact2 = tgamma(5./2. + pars.poly_const)*pars.star_mass/(tgamma(1.+pars.poly_const)*pow(pars.star_radius,3));
        pars.lambda = fact1*pow(fact2,1./pars.poly_const)/pow(pars.star_radius,2);
    }

    //System evolution parameters
    pars.total_time = 10;
    pars.nsteps = 2000;
    pars.dt = pars.total_time/pars.nsteps;

    if(pars.mode == Bidimensional){
        std::vector<Vec2> r;
        std::vector<Vec2> v;
        std::vector<double> m;
        Vec2 pos = Vec2();
        Vec2 vel = Vec2();
        for (int ix = 0; ix<pars.nside; ++ix){
            pos.x = -pars.L/2 + ix*pars.dL;
            for (int iy = 0; iy<pars.nside; ++iy){
                pos.y = -pars.L/2 + iy*pars.dL;
                r.push_back(pos);
                v.push_back(vel);
                m.push_back(pars.particle_mass);
            }
        }
        EvolveIt<Vec2>(pars,&r,&v,&m);
    } else if (pars.mode == Tridimensional){
        std::vector<Vec3> r;
        std::vector<Vec3> v;
        std::vector<double> m;
        Vec3 pos = Vec3();
        Vec3 vel = Vec3();
        for (int ix = 0; ix<pars.nside; ++ix){
            pos.x = -pars.L/2 + ix*pars.dL;
            for (int iy = 0; iy<pars.nside; ++iy){
                pos.y = -pars.L/2 + iy*pars.dL;
                for (int iz = 0; iy<pars.nside; ++iz){
                    pos.z = -pars.L/2 + iz*pars.dL;
                    r.push_back(pos);
                    v.push_back(vel);
                    m.push_back(pars.particle_mass);
                }
            }
        }
        EvolveIt<Vec3>(pars,&r,&v,&m);
    }


}