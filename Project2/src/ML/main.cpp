#include <iostream>
#include "boltzmannmachine.h"
#include "Eigen/Dense"
#include "system.h"
#include "Parameters/parameters.h"
#include "simulation.h"

using namespace std;

int main(int argc, char *argv[])
{
    //Enables Eigen to do matrix operations in parallel
    Eigen::initParallel();

    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    Eigen::ArrayXd test_parameters(8);


    test_parameters << .1,.0,.31,.10,.76,.69,.1,.3;

    System * system = new System();
    Simulation * simulation = new Simulation(system);

    test_parameters = simulation->stochastic_descent(test_parameters);


    std::cout << test_parameters << std::endl;
    delete simulation;
    delete system;

    return 0;
}
