#include "simulation.h"

Simulation::Simulation(int m_MCsteps,int m_size,int m_dim, int m_alpha_min,int m_alpha_max)
{
    trial_function = TrialFunction();
}

void Simulation::initiate(){

}

void Simulation::run(){

}

Particle Simulation::setup_grid(int size, int dim, double a,double m_phi(vec3,double,double)){
    Particle* particle_array[size];
    vec3 r;

    if(dim == 1){
        r = vec3(1,0,0);
    }
    else if(dim == 2){
        r = vec3(1,1,0);
    }
    else if(dim == 3){
        r = vec3(1,1,1);
    }
    else{
        std::cout << "Unknown Number of Dimensions" << std::endl;
    }

    for(int i = 0; i<size; i++){
        particle_array[i] = &Particle(r,a,1,m_phi);
    }
}


