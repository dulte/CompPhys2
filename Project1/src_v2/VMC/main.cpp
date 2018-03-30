#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "Parameters/parameters.h"
#include "simulation.h"
#include "system.h"
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{
    double StartTime = MPI_Wtime();

    Eigen::initParallel();
    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    cout << "MC steps, particles: " << Parameters::MC_cycles << ", " << Parameters::N << endl;

    //Initiates system and simulation
    System * system = new System();
    Simulation * simulation = new Simulation(system);

    simulation->initiate();
    double optimal_alpha = 0.5;//simulation->conjugate_gradient(0.3, 1.);
    //std::cout << "Correct a: " << optimal_alpha <<std::endl;
    //simulation->run();
    //simulation->data_for_derivated();
    simulation->oneBodyDensity(optimal_alpha,0.03,0.,6.);


    double EndTime = MPI_Wtime();
    double TotalTime = EndTime-StartTime;

    cout << "Time = " << TotalTime << endl;

    return 0;
}
