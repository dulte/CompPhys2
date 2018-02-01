#ifndef TRIALFUNCTION_H
#define TRIALFUNCTION_H

#include "Vec3/vec3.h"
#include <math.h>
#include <vector>
#include "particle.h"
#include "potential.h"
#include "Potentials/harmonicoscillator.h"
#include <memory>


class TrialFunction
{
public:
    TrialFunction(std::shared_ptr<Potential>);
    void calculate_trial(std::vector<Particle> p,int,double);
    void calculate_probability();
    void calculate_local_energy(int,int);

    double phi(std::vector<double>, double);


    double get_probability(std::vector<Particle> p,int,double);
    double get_probability_ratio(std::vector<Particle> p, int, int move, double alpha);
    double get_local_energy(int,int);

    std::shared_ptr<Potential> potential;

    double f(std::vector<Particle> p);
    double f_id(std::vector<Particle> p);

private:
    double function_value;
    double function_value_next_step;
    double function_probability;
    double function_probability_next_step;
    double local_energy;
    double (TrialFunction::*f_func)(std::vector<Particle>);
    double beta;
    double a;


};

#endif // TRIALFUNCTION_H
