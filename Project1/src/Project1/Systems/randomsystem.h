#ifndef RANDOMSYSTEM_H
#define RANDOMSYSTEM_H

#include "../system.h"

class RandomSystem : public System
{
public:
    RandomSystem();

    void grid_setup(int,double);
    void propose_step();
    double check_acceptance_and_return_energy();

private:
    double step_size = 0.5;
    double beta = 1;
};

#endif // RANDOMSYSTEM_H
