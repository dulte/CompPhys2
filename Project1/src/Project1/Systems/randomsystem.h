#ifndef RANDOMSYSTEM_H
#define RANDOMSYSTEM_H

#include "../system.h"

class RandomSystem : public System
{
public:
    RandomSystem();
    void grid_setup(int,double);

    double beta = 1;

};

#endif // RANDOMSYSTEM_H
