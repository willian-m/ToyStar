#include <vector>
#include <math.h> //Gamma function
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <sstream>

#include "eos_polytropic.h"
#include "vec3.h"
#include "particle_system.h"
#include "integrator_rk4.h"

#ifdef ROOT_FOUND
#include <typeinfo>
#include <chrono>
#include <thread>
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
//#include "TGraph2D.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#endif

enum HydroType {
    Monodimensional = 1,
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

    bool dynamical_plots_enabled;
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
                <<sys->get_particle(ipart)->get_density()<<std::endl;
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
                <<sys->get_particle(ipart)->get_density()<<std::endl;
        }
    }
    out_file << std::flush;
    out_file.close();
}

void PrintParticleSystem(ParticleSystem<Vec1>* sys,std::string filename){
    std::ofstream out_file(filename,std::ios::trunc);

    //Dumps the data
    if (out_file.is_open() && out_file.good()){
        out_file << "# mass\tx\ty\tz\tvx\tvy\tvz\tdensity"<<std::endl;
        int nparticles = sys->get_nparticles();
        for (int ipart=0; ipart<nparticles;++ipart){
            Particle<Vec1>* p = sys->get_particle(ipart);
            Vec1 pos = p->get_position();
            Vec1 vel = p->get_velocity();
            double x = pos.x; double y = pos.y; double z = .0;
            out_file <<p->get_mass()<<"\t"
                <<x<<"\t"<<y<<"\t"<<z<<"\t"
                <<vel.x<<"\t"<<0<<"\t"<<.0<<"\t"
                <<sys->get_particle(ipart)->get_density()<<std::endl;
        }
    }
    out_file << std::flush;
    out_file.close();
}

void PrintDensityAlongX(ParticleSystem<Vec1>* sys,std::string filename){

    //In addition, gets the density at (z=0, y=0)
    std::ofstream out_density(filename,std::ios::trunc);

    if (out_density.is_open() && out_density.good()){
        out_density << "# r\tdensity"<<std::endl;
        double const max_r=2.;
        double const dr=1.e-2;
        int const nr = max_r/dr;
        for (int ir=0; ir<nr;++ir){
            double r = ir*dr;
            out_density <<r<<"\t"<<sys->get_point_density(Vec1(r))<<std::endl;
        }
    }
    out_density << std::flush;
    out_density.close();
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
            out_density <<r<<"\t"<<sys->get_point_density(Vec2(r,0))<<std::endl;
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
            out_density <<r<<"\t"<<sys->get_point_density(Vec3(r,0,0))<<std::endl;
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
    const clock_t begin_time = clock();

    //Prepare canvas for drawing
    #ifdef ROOT_FOUND

    //Allow plots to be shown
    int argc = 1;
    char* argv[1] = {"ToyStar.exe"};
    TApplication app("Toy Star", &argc, argv);
    
    //Setup analytical solution
    TF1* analytical_sol = new TF1("analytical_sol","[0]*([2]*[2]-x*x)/(4*[1])",0,.75);
    analytical_sol->SetParameter(0, pars.lambda);
    analytical_sol->SetParameter(1, pars.pressure_const);
    analytical_sol->SetParameter(2, pars.star_radius);
    analytical_sol->SetLineColor(kRed);
    
    //Setup axis
    analytical_sol->GetYaxis()->SetRangeUser(0,3);
    analytical_sol->GetXaxis()->SetRangeUser(0,1);
    analytical_sol->GetXaxis()->SetTitle("R");
    analytical_sol->GetYaxis()->SetTitle("#rho");
    
    //Setup plot of particle positions
    TCanvas* cParticlePos = new TCanvas("cParticlePos","Particle Positions", 800, 800);
    TView* view_particles = TView::CreateView(1);
    view_particles->SetRange(-3,-3,-3,3,3,3);
    TPolyMarker3D* particle_pos = nullptr;
    //TGraph2D* particle_pos = nullptr;

    TCanvas* cDensity = new TCanvas("cDensity","Density", 800,800);
    analytical_sol->Draw();
    TGraph* density_graph = nullptr;
    #endif

    for(int istep=0; istep<pars.nsteps;++istep){

        //Evolve the system
        integrator.do_step();
        //Iterate over particles, getting their position, velocity, mass and density
        //Dump it in a text file

        //Do ROOT drawing :)
        #ifdef ROOT_FOUND
        if (pars.dynamical_plots_enabled){
            if (particle_pos != nullptr) delete particle_pos;
            if (density_graph != nullptr) delete density_graph;


            int nparticles = current.get_nparticles();
            cParticlePos->cd();
            particle_pos = new TPolyMarker3D(nparticles);
            //particle_pos = new TGraph2D(nparticles);

            std::vector<double> r, rho;
            for (int ipart=0; ipart<nparticles;++ipart){
                Particle<T>* p = current.get_particle(ipart);
                T pos = p->get_position();
                double x, y, z;
                x = pos.x; y = pos.y;
                z = pos.z;
                particle_pos->SetPoint(ipart,x,y,z);
                r.push_back(sqrt(pow(x,2) + pow(y,2) +pow(z,2) ));

                rho.push_back( p->get_density() );
                //rho.push_back(analytical_sol->Eval(r[ipart])/p->get_density());
            
            }
            
            particle_pos->SetMarkerSize(1);
            particle_pos->SetMarkerColor(kBlue);
            particle_pos->SetMarkerStyle(20);

            particle_pos->Draw("P");
            cParticlePos->Modified(); cParticlePos->Update();

            cDensity->cd();
            density_graph = new TGraph(nparticles,r.data(),rho.data());

            density_graph->SetMarkerStyle(20);
            density_graph->SetMarkerColor(kBlue);

            density_graph->Draw("p");
            cDensity->Modified(); cDensity->Update();

            //app.Run();

            //std::this_thread::sleep_for(std::chrono::milliseconds(200)); //Performance limiter so people can see what is happening :)
        }    
        #endif

        //Progress bar: https://stackoverflow.com/a/14539953/2754579
        float progress = (double) istep/(double) pars.nsteps;
        float time_from_start = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
        float avrg_perf = progress/time_from_start;
        float time_ramaining = (1.-progress)/avrg_perf;
        std::cout << "[";
        int bpos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < bpos) std::cout << "=";
            else if (i == bpos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " % ("<<istep<<"/"<<pars.nsteps<<"). Remaining: "<< time_ramaining/60<<" min\r";
        std::cout.flush();
        if (istep%10 != 0) continue;
        //Create filename
        std::stringstream filename_stream;
        filename_stream << "output_iter_" << istep << ".txt";
        std::string filename;
        filename_stream >> filename;

        //Print output
        PrintParticleSystem(&current, filename);

    }
    PrintDensityAlongX(&current,"density.txt");
}

int main(int argc, char* argv[]){

    //Choose run mode
    ToyStarPars pars;
    pars.mode = Bidimensional;
    //pars.mode = Tridimensional; 
    //pars.mode = Monodimensional;
    pars.dynamical_plots_enabled = true;
    
    //System parameters - initial positions
    pars.L = 2;
    pars.nside = 20*20;//*5*4;
    pars.dL = 2*pars.L/pars.nside;

    //System parameters
    pars.star_mass = 2.;
    pars.star_radius = .75;
    pars.particle_mass = pars.star_mass/pow(pars.nside,(int) pars.mode);
    //const double smoothing = 0.04/sqrt(pow(nside,3)/1000.);
    pars.smoothing = pars.dL;//*1.8;//1.5*0.04/sqrt(pow(nside,3)/1000.);
    //pars.smoothing = .5;
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
    } else if (pars.mode == Monodimensional){
        const double n = pars.poly_const;
        const double fact1 = tgamma(n+1.5)/tgamma(n+1.)/pars.star_radius;
        const double fact2 = pars.star_mass/sqrt(M_PI);
        pars.lambda = pow(fact1*fact2,1./n)*2.*(1.+n)*pars.pressure_const/pow(pars.star_radius,2);
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
                for (int iz = 0; iz<pars.nside; ++iz){
                    pos.z = -pars.L/2 + iz*pars.dL;
                    r.push_back(pos);
                    v.push_back(vel);
                    m.push_back(pars.particle_mass);
                }
            }
        }
        EvolveIt<Vec3>(pars,&r,&v,&m);
    } else if (pars.mode == Monodimensional ){
        std::vector<Vec1> r;
        std::vector<Vec1> v;
        std::vector<double> m;
        Vec1 pos = Vec1();
        Vec1 vel = Vec1();
        for (int ix = 0; ix<pars.nside; ++ix){
            pos.x = -pars.L/2 + ix*pars.dL;
            r.push_back(pos);
            v.push_back(vel);
            m.push_back(pars.particle_mass);
        }
        EvolveIt<Vec1>(pars,&r,&v,&m);
    }


}
