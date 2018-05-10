#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <cmath>
#include "system.h"
#include "Parameters/parameters.h"
#include "DataDump/datadump.h"
#include "Eigen/Dense"


class Simulation
{
public:
    Simulation(System *m_system);
    void initiate();
    void run(int rank,Eigen::ArrayXd &x);


    Eigen::ArrayXd stochastic_descent(Eigen::ArrayXd x_0);
    void calculate_gradient(Eigen::ArrayXd & x,Eigen::ArrayXd &gradient);
    double compute_local_energy_derivative(double alpha);
    void oneBodyDensity(double optimal_alpha, double r_step, double r_min, double r_max);
private:
    System *system;
    double alpha_step;
    double alpha_min;
    double alpha_max;

    double total_energy;

    int MC_cycles;
    int N;

    double energy;


};

#endif // SIMULATION_H
