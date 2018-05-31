#include <iostream>
#include "boltzmannmachine.h"
#include "Eigen/Dense"
#include "system.h"
#include "Parameters/parameters.h"
#include "simulation.h"

using namespace std;

void distribute_weights_and_biases(Eigen::ArrayXd &);

int main(int argc, char *argv[])
{
    //Enables Eigen to do matrix operations in parallel
    Eigen::initParallel();

    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    /*Eigen::ArrayXd test_parameters = Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
    //distribute_weights_and_biases(test_parameters);
    std::cout << "---------------" << std::endl;

    System * system = new System();
    Simulation * simulation = new Simulation(system);

    Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);

    simulation->run(0,done);
    //simulation->run(0,done);

    delete simulation;
    delete system;

    */

    double Ns[5] = {1,2,3,4,5};
    double rates[5] = {0.5,0.1,0.05,0.01,0.005};


    for(int i = 0;i<5;i++){
        for(int j = 0; j<5;j++){
            std::cout << Ns[j] << " "<< rates[i] << std::endl;
            Parameters::N = Ns[j];
            Parameters::learning_rate = rates[i];
            System * system = new System();
            Simulation * simulation = new Simulation(system);
            Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
            distribute_weights_and_biases(test_parameters);
            std::cout << "ehi " << Parameters::N << std::endl;
            Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
            delete simulation;
            delete system;
        }

    }

    /*
    double sigma = 0.5;
    for(int i = 0;i<50;i++){
        Parameters::sigma = sigma;
        System * system = new System();
        Simulation * simulation = new Simulation(system);
        Eigen::ArrayXd test_parameters = Eigen::ArrayXd::Zero(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
        distribute_weights_and_biases(test_parameters);
        std::cout << "ehi " << Parameters::N << std::endl;
        Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
        simulation->run(0,done);
        simulation->run(0,done);
        delete simulation;
        delete system;

        sigma += 0.02;
    }
    */

    /*
    //double dx[7] = {1.5,1.25,1,0.75,0.5,0.25,0.1};
    double dx[7] = {1,0.5,0.1,0.05,0.01,0.005,0.001};
    Eigen::ArrayXd test_parameters = Eigen::ArrayXd::Zero(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
    distribute_weights_and_biases(test_parameters);
    //Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
    for(int i = 0;i<7;i++){
        Parameters::dx = dx[i];
        System * system = new System();
        Simulation * simulation = new Simulation(system);
        Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
        simulation->run(0,done);
        delete simulation;
        delete system;
    }
    */

    return 0;
}

void distribute_weights_and_biases(Eigen::ArrayXd & array){
    int size = array.size();


    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<double> distribution(0,0.5);
    for(int i = 0;i<size;i++){
        array[i] = distribution(gen);
    }
}
