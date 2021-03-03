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

int main(int argc, char* argv[]){

    //Choose run mode
    HydroType mode = Bidimensional;

    //System parameters - initial positions
    const double L = 2*.75;
    const int nside = 20;
    const double dL = L/nside;

    //System parameters
    const double star_mass = 2.;
    const double star_radius = .75;
    const double particle_mass = star_mass/pow(nside,(int) mode);
    //const double smoothing = 0.04/sqrt(pow(nside,3)/1000.);
    const double smoothing = dL;//*1.8;//1.5*0.04/sqrt(pow(nside,3)/1000.);
    const double damping = 1.;

    //Allocates eos
    const double pressure_const = 1;
    const double poly_const = 1;
    EOSPolytropic* eos = new EOSPolytropic(pressure_const,poly_const);

    //Get attractive constant
    double lambda = 20;
    if (mode == Bidimensional){
        const double fact1 = 2*pressure_const*pow(M_PI,-1./poly_const);
        const double fact2 = pow(star_mass*(1+poly_const)/pow(star_radius,2),1.+1./poly_const);
        lambda = fact1*fact2/star_mass;
    } else if (mode == Tridimensional){
        const double fact1 = 2*pressure_const*(1+poly_const)*pow(M_PI,-3./(2*poly_const));
        const double fact2 = tgamma(5./2. + poly_const)*star_mass/(tgamma(1.+poly_const)*pow(star_radius,3));
        lambda = fact1*pow(fact2,1./poly_const)/pow(star_radius,2);
    }

    //System evolution parameters
    const double total_time = 16;
    const int nsteps = 400;
    const double dt = total_time/nsteps;

    if (mode == Tridimensional){

        //Creates particles
        std::vector<Vec3> r;
        std::vector<Vec3> v;
        std::vector<double> m;
        Vec3 pos(.0, .0, .0);
        Vec3 vel(.0, .0, .0);
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

        //Each particle behaves as an damped HO. We use our knowledge of the 
        //anallytic solution to estimate how much time we must wait to get 
        //to the stationary solution 

        //Print parameters
        std::cout << "SPH parameters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Number of particles..: " << pow(nside,3) << std::endl;
        std::cout << "Particle mass........: " << particle_mass << std::endl;
        std::cout << "Smoothing (h)........: " << smoothing << std::endl;

        std::cout << "Star paramneters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Damping..............: " << damping << std::endl;
        std::cout << "Gravity..............: " << lambda << std::endl;
        std::cout << "Star mass............: " << star_mass << std::endl;
        std::cout << "Star radius..........: " << star_radius << std::endl;


        std::cout << "EOS parameters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Pressure constant....: " << pressure_const << std::endl;
        std::cout << "Polytropic constant..: " << poly_const << std::endl;

        std::cout << "Simulation parameters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Simulation time......: " << total_time << std::endl;
        std::cout << "Number of steps......: " << nsteps << std::endl;
        std::cout << "Step size............: " << dt << std::endl;

        //Init particle system
        ParticleSystem<Vec3> current(r,v,m,smoothing,lambda,damping,eos);
        ParticleSystem<Vec3> next(r,v,m,smoothing,lambda,damping,eos);
        ParticleSystem<Vec3> buffer(r,v,m,smoothing,lambda,damping,eos);

        //Init integrator
        IntegratorRK4<Vec3> integrator(&current,&next,&buffer,dt);

        //Simulate it!    
        int barWidth = 70; //Width of progress bar
        for(int istep=0; istep<nsteps;++istep){

            //Evolve the system
            integrator.do_step();
            //Iterate over particles, getting their position, velocity, mass and density
            //Dump it in a text file

            //a) Create filename
            std::stringstream filename_stream;
            filename_stream << "output_iter_" << istep << ".txt";
            std::string filename;
            filename_stream >> filename;

            std::ofstream out_file(filename,std::ios::trunc);

            //Dumps the data
            if (out_file.is_open() && out_file.good()){
                out_file << "# mass\tx\ty\tz\tvx\tvy\tvz\tdensity"<<std::endl;
                int nparticles = current.get_nparticles();
                for (int ipart=0; ipart<nparticles;++ipart){
                    Particle<Vec3>* p = current.get_particle(ipart);
                    Vec3 pos = p->get_position();
                    Vec3 vel = p->get_velocity();
                    double x = pos.x; double y = pos.y; double z = pos.z;
                    out_file <<p->get_mass()<<"\t"
                             <<x<<"\t"<<y<<"\t"<<z<<"\t"
                             <<vel.x<<"\t"<<vel.y<<"\t"<<vel.z<<"\t"
                             <<current.get_density(Vec3(x,y,z))<<std::endl;
                }
            }
            out_file << std::flush;
            out_file.close();


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



        //In addition, gets the density at (z=0, y=0)
        std::ofstream out_density("density.txt",std::ios::trunc);

        if (out_density.is_open() && out_density.good()){
            out_density << "# r\tdensity"<<std::endl;
            double const max_r=2.;
            double const dr=1.e-2;
            int const nr = max_r/dr;
            for (int ir=0; ir<nr;++ir){
                double r = ir*dr;
                out_density <<r<<"\t"<<current.get_density(Vec3(r,0,0))<<std::endl;
            }
        }
        out_density << std::flush;
        out_density.close();
    } else
    if (mode == Bidimensional) {
        

        //Creates particles
        std::vector<Vec2> r;
        std::vector<Vec2> v;
        std::vector<double> m;
        Vec2 pos(.0, .0);
        Vec2 vel(.0, .0);
        for (int ix = 0; ix<nside; ++ix){
            pos.x = -L/2 + ix*dL;
            for (int iy = 0; iy<nside; ++iy){
                pos.y = -L/2 + iy*dL;
                r.push_back(pos);
                v.push_back(vel);
                m.push_back(particle_mass);
                }
            }
        

        //Each particle behaves as an damped HO. We use our knowledge of the 
        //anallytic solution to estimate how much time we must wait to get 
        //to the stationary solution 

        //Print parameters
        std::cout << "SPH parameters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Number of particles..: " << pow(nside,2) << std::endl;
        std::cout << "Particle mass........: " << particle_mass << std::endl;
        std::cout << "Smoothing (h)........: " << smoothing << std::endl;

        std::cout << "Star paramneters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Damping..............: " << damping << std::endl;
        std::cout << "Gravity..............: " << lambda << std::endl;
        std::cout << "Star mass............: " << star_mass << std::endl;
        std::cout << "Star radius..........: " << star_radius << std::endl;


        std::cout << "EOS parameters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Pressure constant....: " << pressure_const << std::endl;
        std::cout << "Polytropic constant..: " << poly_const << std::endl;

        std::cout << "Simulation parameters" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Simulation time......: " << total_time << std::endl;
        std::cout << "Number of steps......: " << nsteps << std::endl;
        std::cout << "Step size............: " << dt << std::endl;

        //Init particle system
        ParticleSystem<Vec2> current(r,v,m,smoothing,lambda,damping,eos);
        ParticleSystem<Vec2> next(r,v,m,smoothing,lambda,damping,eos);
        ParticleSystem<Vec2> buffer(r,v,m,smoothing,lambda,damping,eos);

        //Init integrator
        IntegratorRK4<Vec2> integrator(&current,&next,&buffer,dt);

        //Simulate it!    
        int barWidth = 70; //Width of progress bar
        for(int istep=0; istep<nsteps;++istep){

            //Evolve the system
            integrator.do_step();
            //Iterate over particles, getting their position, velocity, mass and density
            //Dump it in a text file

            //a) Create filename
            std::stringstream filename_stream;
            filename_stream << "output_iter_" << istep << ".txt";
            std::string filename;
            filename_stream >> filename;

            std::ofstream out_file(filename,std::ios::trunc);

            //Dumps the data
            if (out_file.is_open() && out_file.good()){
                out_file << "# mass\tx\ty\tvx\tvy\tdensity"<<std::endl;
                int nparticles = current.get_nparticles();
                for (int ipart=0; ipart<nparticles;++ipart){
                    Particle<Vec2>* p = current.get_particle(ipart);
                    Vec2 pos = p->get_position();
                    Vec2 vel = p->get_velocity();
                    double x = pos.x; double y = pos.y;
                    out_file <<p->get_mass()<<"\t"
                             <<x<<"\t"<<y<<"\t"
                             <<vel.x<<"\t"<<vel.y<<"\t"
                             <<current.get_density(Vec2(x,y))<<std::endl;
                }
            }
            out_file << std::flush;
            out_file.close();


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



        //In addition, gets the density at (z=0, y=0)
        std::ofstream out_density("density.txt",std::ios::trunc);

        if (out_density.is_open() && out_density.good()){
            out_density << "# r\tdensity"<<std::endl;
            double const max_r=2.;
            double const dr=1.e-2;
            int const nr = max_r/dr;
            for (int ir=0; ir<nr;++ir){
                double r = ir*dr;
                out_density <<r<<"\t"<<current.get_density(Vec2(r,0))<<std::endl;
            }
        }
        out_density << std::flush;
        out_density.close();

    }
    

    




}