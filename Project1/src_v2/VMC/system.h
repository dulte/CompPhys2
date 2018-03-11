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
    Eigen::MatrixXd quantum_force_matrix;
    Eigen::MatrixXd quantum_force_matrix_new;

    const Eigen::MatrixXd get_position() const;

    int acceptance = 0;
    void make_grid(double m_alpha);
    void make_move_and_update_non_interacting(const int move);
    void make_move_and_update_interacting(const int move);
    void make_move_and_update(const int move);
    void update();
    double check_acceptance_and_return_energy(int);


    void update_wavefunction(const int move);
    double calculate_energy_noninteracting();
    double calculate_energy();
    double calculate_energy_interacting();
    double greens_function_ratio(int move);
    void update_next_distance(int);

    double expectation_local_energy;
    double expectation_derivative;
    double expectation_derivative_energy;
    double expectation_local_energy_squared;

    double local_energy_interacting();
    double udiv(int,int);
    double udivdiv(int,int);


    double update_wavefunction_interacting(int);
    double f(double);
private:
    //Saves all the variables from the parameters to save time
    const int N = Parameters::N;
    const int dimension = Parameters::dimension;
    const double beta = Parameters::beta;
    const double dx = Parameters::dx;
    const double a = Parameters::a;
    const double omega = Parameters::omega;
    const double omega_z = Parameters::omega_z;
    const double D = Parameters::D;
    double h;
    double wavefunction_value_plus;
    double wavefunction_value_minus;


    void (System::*make_move)(const int);
    double (System::*compute_energy)();


    //Vector and double used for holding temp values
    Eigen::VectorXd temp_r;
    Eigen::VectorXd temp_r2;
    Eigen::VectorXd temp_N_vec;
    //double temp_value;
    //double temp_value2;



    //Functions for calulations involving wavefunction
    double alpha;
    double wavefunction_value;
    double wavefunction_probability;
    double local_energy;
    double kinetic_energy;
    double potential_energy;

    double phi_exponant(const Eigen::VectorXd &r);
    double get_wavefunction();
    double get_probability_ratio(int move);
    double get_probability();
    double get_local_energy();
    void update_probability_ratio();

    void quantum_force(int);

    //Functions for the use of normal distributions
    std::random_device rd;
    std::mt19937_64 gen;
    std::normal_distribution<double> distribution;


};

#endif // SYSTEM_H
