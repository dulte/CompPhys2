#include "randomsystem.h"

RandomSystem::RandomSystem()
{
    dimension = 1;
    trial_function = TrialFunction();
}

void RandomSystem::grid_setup(int size, double start_alpha)
{
    vec3 r;
    r = vec3(1,0,0);

    for(int i = 0; i<size; i++){
        phi_values.push_back(trial_function.phi(r,start_alpha,beta));
        particles.push_back(Particle(r,start_alpha,1));
    }

}

