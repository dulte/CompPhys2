#ifndef SIMULATION_H
#define SIMULATION_H
#include "system.h"
#include "Parameters/parameters.h"
#include "DataDump/datadump.h"


class Simulation
{
public:
    Simulation(System *m_system);
    void initiate();
    void run();

private:
    System *system;
    double alpha_step;
    double alpha_min;
    double alpha_max;

    int MC_cycles;
    int N;

    double energy;


};

#endif // SIMULATION_H
