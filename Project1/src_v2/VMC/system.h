#ifndef SYSTEM_H
#define SYSTEM_H
#include "Eigen/Dense"
#include "Parameters/parameters.h"
#include <random>
#include <time.h>
#include <math.h>
#include <iostream>

class System
{
public:
    System();

    Eigen::MatrixXd r;
    Eigen::MatrixXd next_r;

    Eigen::MatrixXd distance;
    Eigen::MatrixXd next_distance;
    Eigen::VectorXd quantum_force_vector;
    Eigen::VectorXd quantum_force_vector_new;

    const Eigen::MatrixXd get_position() const;

    void make_grid(double m_alpha);
    void make_move_and_update(int move);
    void update();
    double check_acceptance_and_return_energy(int);

    double expectation_local_energy;
    double expectation_derivative;
    double expectation_derivative_energy;
    double expectation_local_energy_squared;
    void update_expectation();

    void update_next_distance(int move);

    void quantum_force(int move);
    double greens_function_ratio(int move);
    void distribute_particles_interacting();
    void distribute_particles_noninteracting();
private:
    //Saves all the variables from the parameters to save time
    const int N = Parameters::N;
    const int dimension = Parameters::dimension;
    const double beta = Parameters::beta;
    const double dx = Parameters::dx;
    const double omega = Parameters::omega;
    const double omega_z = Parameters::omega_z;
    const double a = Parameters::a;
    const double D = Parameters::D;


    //Vector and double used for holding temp values
    Eigen::VectorXd temp_r;
    double temp_value;
    //double temp_value2;
    Eigen::MatrixXd *r_temp_pointer;



    void (System::*wavefunction_function_pointer)(const int);
    double (System::*compute_energy_numeric)();
    double (System::*compute_local_energy)();


    //Functions for calulations involving wavefunction
    double alpha;
    double wavefunction_value;
    double wavefunction_probability;
    double local_energy;

    double phi_exponant(const Eigen::VectorXd &r);
    double get_wavefunction();
    double get_probability_ratio(int move);
    double get_probability();
    double get_local_energy();

    void update_wavefunction(const int move);
    void update_probability_ratio();

    double udiv(int,int);
    double udivdiv(int,int);
    double get_local_energy_interacting();
    double update_wavefunction_interacting_f(const int);
    void update_wavefunction_interacting(const int);
    void update_wavefunction_noninteracting(const int move);
    double get_local_energy_noninteracting();
    double f(double);



    std::random_device rd;
    std::mt19937_64 gen;
    std::normal_distribution<double> distribution;




};

#endif // SYSTEM_H
