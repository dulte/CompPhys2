#ifndef SIMULATION_H
#define SIMULATION_H

#include "particle.h"
#include "trialfunction.h"
#include "system.h"
#include "potential.h"
#include "Systems/randomsystem.h"
#include "DataDump/datadump.h"
#include <iostream>
#include <vector>
#include <memory>


class Simulation
{
public:
    Simulation(System *m_system);
    void initiate();
    void run(int);

private:
    int MCsteps;
    double alpha_min;
    double alpha_max;
    double alpha_step;
    int size;
    System *system;
    Potential *potential;
    TrialFunction *trial_function;

    //DataDump<double> dump;
};

#endif // SIMULATION_H
