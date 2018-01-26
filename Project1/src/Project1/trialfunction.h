#ifndef TRIALFUNCTION_H
#define TRIALFUNCTION_H

#include "Vec3/vec3.h"
#include <math.h>
#include "particle.h"


class TrialFunction
{
public:
    TrialFunction();
    void calculate_trial(Particle *p,int,double,double);
    void calculate_probability();
    void calculate_local_energy(int,int);

    static double phi(vec3,double,double);

    double get_probability(Particle *p,int,double,double);
    double get_local_energy(int,int);

private:
    double function_value;
    double function_probability;
    double local_energy;
};

#endif // TRIALFUNCTION_H
