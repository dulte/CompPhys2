#include <iostream>
#include "boltzmannmachine.h"
#include "Eigen/Dense"
#include "system.h"
#include "Parameters/parameters.h"

using namespace std;

int main(int argc, char *argv[])
{
    //Enables Eigen to do matrix operations in parallel
    Eigen::initParallel();

    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    Eigen::ArrayXd test_parameters(8);

    test_parameters << 1,2,3,4,5,6,7,8;

    System syst;
    syst.make_grid(test_parameters);

    return 0;
}
