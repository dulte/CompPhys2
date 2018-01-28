#include "randomsystem.h"
#include <iostream>

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
    vec3 r;
    r = vec3(1,0,0);


    for(int i = 0; i<size; i++){
        phi_values.push_back(trial_function->phi(r,start_alpha,beta));
        particles.push_back(Particle(r,start_alpha,1));
    }

}

void RandomSystem::propose_step(){
    vec3 step;

    for(int i = 0; i< size; i++){
        step = vec3((float)rand()/RAND_MAX,(float)rand()/RAND_MAX,(float)rand()/RAND_MAX);
        particles[i].next_r = particles[i].r + step_size*step;
    }
}

double RandomSystem::check_acceptance_and_return_energy(){

    double delta_energi = 0;
    double r = (float)rand()/RAND_MAX;
    double acceptance_probability = 0;

    std::cout << "hei" << std::endl;

    trial_function->get_probability(particles,size,alpha,beta);

    for(int i = 0; i< size; i++){
        r = (float)rand()/RAND_MAX;
        acceptance_probability = trial_function->get_probability_ratio(particles,size,i,alpha,beta);
        particles[i].r.print();
        if(acceptance_probability >= r){
            particles[i].accept_step();
            trial_function->get_probability(particles,size,alpha,beta); //Not sure if should be here
        }
    }
}

