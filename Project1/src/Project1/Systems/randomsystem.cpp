#include "randomsystem.h"
#include <iostream>
#include <math.h>

RandomSystem::RandomSystem(std::shared_ptr<TrialFunction> m_trial)
{
    dimension = 1;
    trial_function = m_trial;//TrialFunction();

    std::cout << "MC_cylces: " << Parameters::MC_cycles << std::endl;
}

void RandomSystem::update_alpha(double m_alpha)
{
    alpha = m_alpha;
}

void RandomSystem::grid_setup(int m_size, double start_alpha)
{
    alpha = start_alpha;
    size = Parameters::N;
    std::vector<double> r;
    for(int i = 0; i<dimension; i++){
        r.push_back(0.5);
    }



    for(int i = 0; i<size; i++){
        phi_values.push_back(trial_function->phi(r,start_alpha));
        particles.push_back(Particle(r,start_alpha,1));
    }

}

void RandomSystem::propose_step(){


    for(int i = 0; i< size; i++){
        for(int j = 0; j<dimension; j++){
            particles[i].next_r[j] = particles[i].r[j] + step_size*((float)rand()/RAND_MAX - 0.5);
        }
    }
}

double RandomSystem::check_acceptance_and_return_energy(){

    double delta_energy = 0;
    double r = (float)rand()/RAND_MAX;
    double acceptance_probability = 0;


    trial_function->get_probability(particles,size,alpha);

    for(int i = 0; i< size; i++){
        r = (float)rand()/RAND_MAX;
        acceptance_probability = trial_function->get_probability_ratio(particles,size,i,alpha);
        if(acceptance_probability >= r){
            particles[i].accept_step();
            std::cout << "Pos: " << particles[i].r[0] << std::endl;
            trial_function->get_probability(particles,size,alpha); //Not sure if should be here
        }

        delta_energy = delta_energy + trial_function->get_local_energy(1,1);
    }
    return delta_energy;
}

