#ifndef SYSTEM_H
#define SYSTEM_H

#include "particle.h"
#include "trialfunction.h"
#include "Vec3/vec3.h"
#include <iostream>
#include <list>

class System
{
public:
    //System();
    ~System();
    virtual void grid_setup(int,double){};

    std::list<Particle> particles;
    std::list<double> phi_values;
    TrialFunction trial_function;
    int dimension = 3;
};

#endif // SYSTEM_H
