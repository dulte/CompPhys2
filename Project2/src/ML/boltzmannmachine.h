#pragma once
#include "Eigen/Dense"
#include <cmath>

class BoltzmannMachine
{
public:
    BoltzmannMachine(int visibleNodes,int hiddenNode);

    Eigen::VectorXd feedForward(Eigen::VectorXd input);
    Eigen::VectorXd feedBackward(Eigen::VectorXd input);

private:
    Eigen::VectorXd visibleBias;
    Eigen::VectorXd hiddenBias;

    Eigen::MatrixXd weights;

    Eigen::VectorXd activation(const Eigen::VectorXd &);
};

