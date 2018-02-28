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

    Eigen::MatrixXd *r;
    Eigen::MatrixXd *next_r;

    Eigen::MatrixXd distance;
    Eigen::MatrixXd next_distance;

    const Eigen::MatrixXd get_position() const;

    void make_grid(double m_alpha);
    void make_move_and_update(int move);
    void update();
    double check_acceptance_and_return_energy();

private:
    //Saves all the variables from the parameters to save time
    const int N = Parameters::N;
    const int dimension = Parameters::dimension;
    const double beta = Parameters::beta;
    const double dx = Parameters::dx;

    //Vector and double used for holding temp values
    //Eigen::VectorXd temp_r;
    //double temp_value;
    //double temp_value2;
    Eigen::MatrixXd *r_temp_pointer;


    //Functions for calulations involving wavefunction
    double alpha;
    double wavefunction_value;
    double wavefunction_probability;
    double local_energy;

    double phi_exponant(const Eigen::VectorXd &r);
    double get_wavefunction();
    double get_probability_ratio();
    double get_probability();
    double get_local_energy();

    void update_wavefunction();
    void update_probability_ratio();




};

#endif // SYSTEM_H
