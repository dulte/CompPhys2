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
    size = m_size;
    std::vector<double> r;
    r.push_back(0.5);


    for(int i = 0; i<size; i++){
        phi_values.push_back(trial_function->phi(r,start_alpha,beta));
        particles.push_back(Particle(r,start_alpha,1));
    }

}

void RandomSystem::propose_step(){


    for(int i = 0; i< size; i++){
        for(int j = 0; j<dimension; j++){
            particles[i].next_r[j] = particles[i].r[j] + step_size*(2*(float)rand()/RAND_MAX - 1);
        }
    }
}

double RandomSystem::check_acceptance_and_return_energy(){

    double delta_energi = 0;
    double r = (float)rand()/RAND_MAX;
    double acceptance_probability = 0;


    trial_function->get_probability(particles,size,alpha,beta);

    for(int i = 0; i< size; i++){
        r = (float)rand()/RAND_MAX;
        acceptance_probability = trial_function->get_probability_ratio(particles,size,i,alpha,beta); 
        if(acceptance_probability >= r){
            particles[i].accept_step();
            std::cout << "Pos: " << particles[i].r[0] << std::endl;
            trial_function->get_probability(particles,size,alpha,beta); //Not sure if should be here
        }
        //Remember to update delta_E

        delta_energi = delta_energi + trial_function->get_local_energy(1,1);
    }
    return delta_energi;
}

