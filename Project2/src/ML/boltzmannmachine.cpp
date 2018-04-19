#include "boltzmannmachine.h"
#include <iostream>

BoltzmannMachine::BoltzmannMachine(int visibleNodes, int hiddenNodes)
{
    visibleBias = Eigen::VectorXd::Random(visibleNodes);
    hiddenBias = Eigen::VectorXd::Random(hiddenNodes);
    weights =  Eigen::MatrixXd::Random(hiddenNodes,visibleNodes);

    /*weights << 1,2,
            3,4,
            5,6,
            7,8;*/

}

Eigen::VectorXd BoltzmannMachine::activation(const Eigen::VectorXd & input){
    Eigen::ArrayXd x = (Eigen::ArrayXd) input;
    return (Eigen::VectorXd) Eigen::inverse(1+(-x).exp());
}

Eigen::VectorXd BoltzmannMachine::feedForward(Eigen::VectorXd input){
    return activation(weights.transpose()*input);
}

Eigen::VectorXd BoltzmannMachine::feedBackward(Eigen::VectorXd hidden){
    return activation(weights*hidden);
}


