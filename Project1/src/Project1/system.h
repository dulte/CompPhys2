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
#include <memory>

using namespace std;

class System
{
public:
    System();


    virtual void update_alpha(double){};

    virtual void grid_setup(int,double){cout << "Nei!" << endl;};
    virtual void propose_step(){};
    virtual double check_acceptance_and_return_energy(){};

    std::vector<Particle> particles;

    std::vector<std::vector<double>> particle_r;
    std::vector<std::vector<double>> particle_r_next;
    std::vector<double> phi_values;
    TrialFunction *trial_function;
    int dimension = Parameters::dimension;
    double alpha;
    int size;

    std::vector<double> rs;
    std::vector<double> get_postions();
};

#endif // SYSTEM_H
