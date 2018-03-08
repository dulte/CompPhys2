#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Parameters/parameters.h"
#include "simulation.h"
#include "system.h"

using namespace std;

int main(int argc, char *argv[])
{

    Eigen::initParallel();
    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    cout << "MC steps, particles: " << Parameters::MC_cycles << ", " << Parameters::N << endl;

    //Initiates system and simulation
    System * system = new System();
    Simulation * simulation = new Simulation(system);

    simulation->initiate();
    //simulation->run();
    std::cout << "Correct a: " << simulation->conjugate_gradient(0.5, 1e-7)<<std::endl;

    delete simulation;
    delete system;

    return 0;
}
