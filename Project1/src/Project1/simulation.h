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
    Simulation(std::shared_ptr<System> m_system);
    void initiate();
    void run(int);

private:
    int MCsteps;
    double alpha_min;
    double alpha_max;
    double alpha_step;
    int size;
    std::shared_ptr<System> system;
    std::shared_ptr<Potential> potential;
    std::shared_ptr<TrialFunction> trial_function;

    //DataDump<double> dump;
};

#endif // SIMULATION_H
