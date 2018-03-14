#ifndef SIMULATION_H
#define SIMULATION_H
#include "system.h"
#include "Parameters/parameters.h"
#include "DataDump/datadump.h"
#include <string>


class Simulation
{
public:
    Simulation(System *m_system);
    void initiate();
    void run(int);

    double compute_local_energy_derivative(double alpha);
    double conjugate_gradient(double alpha_0, double b);
    void check_derivative_of_energy();
private:
    System *system;
    double alpha_step;
    double alpha_min;
    double alpha_max;

    int MC_cycles;
    int N;

    double total_energy;
    double energy;
    double energy_numerical;




};

#endif // SIMULATION_H
