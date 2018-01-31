#include "trialfunction.h"
#include <iostream>

TrialFunction::TrialFunction(std::shared_ptr<Potential> m_potential)
{
    potential = m_potential;
}

void TrialFunction::calculate_trial(std::vector<Particle> p, int size, double alpha, double beta)
{
    double val = 1;
    std::vector<double> r;

    for(int i = 0; i<size;i++){
        r = p[i].r;
        val *= phi(r,alpha,beta);
    }

    function_value = val;
}

void TrialFunction::calculate_probability(){
    function_probability = function_value*function_value;
    //function_probability_next_step = function_value_next_step*function_value_next_step;
}

void TrialFunction::calculate_local_energy(int n,int dim){
    local_energy = dim/2.0*n; // And external!!!
}

double TrialFunction::phi(std::vector<double> r, double alpha, double beta)
{
    int size = r.size();
    double val = 0;
    for(int i = 0;i<size;i++){
        if(i == 2){
            val += beta*r[i]*r[i];
        }
        else{
            val += r[i]*r[i];
        }
    }

    return exp(-alpha*val);
}


double TrialFunction::get_probability(std::vector<Particle> p,int size,double alpha,double beta){
    calculate_trial(p,size,alpha,beta);
    calculate_probability();
    return function_probability;
}

double TrialFunction::get_probability_ratio(std::vector<Particle> p,int size,int move, double alpha, double beta){
    double val = 1;
    std::vector<double> r;

    for(int i = 0; i<size;i++){
        if(i == move){
            r = p[i].next_r;
        }
        else{
            r = p[i].r;
        }

        val *= phi(r,alpha,beta); //Can be optimized since this is a product of exp

    }
    return val*val/function_probability;
}



double TrialFunction::get_local_energy(int n, int dim){
    calculate_local_energy(n,dim);
    return local_energy;
}