#ifndef SIMULATION_H
#define SIMULATION_H

#include "particle.h"
#include "trialfunction.h"
#include <iostream>


class Simulation
{
public:
    Simulation(int,int,int,int,int);
    Particle setup_grid(int,int,double,double m_phi(vec3,double,double));
    void initiate();
    void run();

private:
    int m_MCsteps;
    double* alphas;
    Particle* particles;
    TrialFunction trial_function;
};

#endif // SIMULATION_H
