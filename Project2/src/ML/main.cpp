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

    Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);

    std::cout << test_parameters << std::endl;
    std::cout << "---------------" << std::endl;
    //test_parameters << -.1,.0,-.31,.10,.76,.69,-.1,.3;

    System * system = new System();
    Simulation * simulation = new Simulation(system);

    Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);


    std::cout <<  test_parameters - done << std::endl;
    delete simulation;
    delete system;

    return 0;
}
