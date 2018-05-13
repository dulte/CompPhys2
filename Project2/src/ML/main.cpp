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

    /*Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);

    std::cout << "---------------" << std::endl;

    System * system = new System();
    Simulation * simulation = new Simulation(system);

    Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);

    std::cout <<  test_parameters - done << std::endl;


    delete simulation;
    delete system;*/



    double Ns[5] = {5,4,3,2,1};
    double rates[4] = {0.5,0.1,0.05,0.01};

    for(int i = 0;i<4;i++){
        for(int j = 0; j<5;j++){
            Parameters::N = Ns[j];
            Parameters::learning_rate = rates[i];
            System * system = new System();
            Simulation * simulation = new Simulation(system);
            Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
            std::cout << "ehi " << Parameters::N << std::endl;
            Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
            delete simulation;
            delete system;
        }

    }


    return 0;
}
