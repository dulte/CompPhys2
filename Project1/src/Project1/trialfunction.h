#ifndef TRIALFUNCTION_H
#define TRIALFUNCTION_H

#include "Vec3/vec3.h"
#include <math.h>
#include <vector>
#include "particle.h"


class TrialFunction
{
public:
    TrialFunction();
    void calculate_trial(std::vector<Particle> p,int,double,double);
    void calculate_probability();
    void calculate_local_energy(int,int);

    static double phi(vec3,double,double);

    double get_probability(std::vector<Particle> p,int,double,double);
    double get_probability_ratio(std::vector<Particle> p, int, int move, double alpha, double beta);
    double get_local_energy(int,int);

private:
    double function_value;
    double function_value_next_step;
    double function_probability;
    double function_probability_next_step;
    double local_energy;
};

#endif // TRIALFUNCTION_H
