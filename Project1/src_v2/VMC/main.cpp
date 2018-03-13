#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
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
    simulation->run();

    double EndTime = MPI_Wtime();
        double TotalTime = EndTime-StartTime;

    cout << "Time = " << TotalTime << endl;

    return 0;
}
