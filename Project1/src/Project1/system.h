#ifndef SYSTEM_H
#define SYSTEM_H

#include "particle.h"
#include "trialfunction.h"
#include "Vec3/vec3.h"
#include <iostream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <time.h>

using namespace std;

class System
{
public:
    //System();
    ~System();

    void update_alpha(double);

    virtual void grid_setup(int,double){}
    virtual void propose_step(){}
    virtual double check_acceptance_and_return_energy(){return 0;}

    std::vector<Particle> particles;
    std::vector<double> phi_values;
    TrialFunction trial_function;
    int dimension = 3;
    double alpha;
    int size;
};

#endif // SYSTEM_H
