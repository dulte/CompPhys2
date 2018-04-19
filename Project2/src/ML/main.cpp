#include <iostream>
#include "boltzmannmachine.h"
#include "Eigen/Dense"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;

    BoltzmannMachine BM(4,2);
    Eigen::VectorXd input(2);
    input << .1,.2;

    Eigen::VectorXd hidden(4);
    hidden << .1,.2,.3,.4;
    std::cout << BM.feedForward(input) << std::endl;
    std::cout << BM.feedBackward(hidden) << std::endl;

    return 0;
}
