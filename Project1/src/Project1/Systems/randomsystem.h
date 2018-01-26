#ifndef RANDOMSYSTEM_H
#define RANDOMSYSTEM_H

#include "../system.h"

class RandomSystem : public System
{
public:
    RandomSystem();
    void grid_setup(int,double);

    double beta = 1;

    void propose_step();
    double check_acceptance_and_return_energy();
    double step_size = 0.5;

};

#endif // RANDOMSYSTEM_H
