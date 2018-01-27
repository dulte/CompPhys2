#ifndef SIMULATION_H
#define SIMULATION_H

#include "particle.h"
#include "trialfunction.h"
#include "system.h"
#include "Systems/randomsystem.h"
#include <iostream>
#include <vector>


class Simulation
{
public:
    Simulation();
    void initiate(int,int,int,int);
    void run(int);

private:
    int MCsteps;
    double alpha_min;
    double alpha_max;
    double alpha_step;
    int size;
    RandomSystem system;


};

#endif // SIMULATION_H
