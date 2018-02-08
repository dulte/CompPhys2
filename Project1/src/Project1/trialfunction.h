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
    void quantum_force(std::vector<Particle> p, double alpha);
    double greens_function_ratio(std::vector<Particle> p, double alpha, int chosen_particle);
    double greens_function_ratio_id(std::vector<Particle> p, double alpha, int chosen_particle);
    void quantum_force_new(std::vector<Particle> p, double alpha, int chosen_particle);
    void allocate_empty_arrays();
private:
    double function_value;
    double function_value_next_step;
    double function_probability;
    double function_probability_next_step;
    double local_energy;
    double (TrialFunction::*f_func)(std::vector<Particle>);
    double (TrialFunction::*greens_function_ratio_func)(std::vector<Particle>, double, int);
    double beta;
    double dx;
    double a;
    int dimension;
    int N;
    double D;
    std::vector<std::vector<double>> quantum_force_matrix;
    std::vector<std::vector<double>> quantum_force_matrix_new;
    std::vector<std::vector<double>> distance_matrix;
    std::vector<std::vector<double>> distance_matrix_new;

};

#endif // TRIALFUNCTION_H
