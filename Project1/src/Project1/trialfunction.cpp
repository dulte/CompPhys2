#include "trialfunction.h"

TrialFunction::TrialFunction()
{

}

void TrialFunction::calculate_trial(std::vector<Particle> p, int size, double alpha, double beta)
{
    double val = 0;
    vec3 r;

    for(int i = 0; i<size;i++){
        r = p[i].r;
        val *= phi(r,alpha,beta);
    }

    function_value = val;
}

void TrialFunction::calculate_probability(){
    function_probability = function_value*function_value;
    function_probability_next_step = function_value_next_step*function_value_next_step;
}

void TrialFunction::calculate_local_energy(int n,int dim){
    local_energy = dim/2.0*n;
}

double TrialFunction::phi(vec3 r, double alpha, double beta)
{
    return exp(-alpha*(r[0]*r[0] + r[1]*r[1] + beta*r[2]*r[2]));
}


double TrialFunction::get_probability(std::vector<Particle> p,int size,double alpha,double beta){
    calculate_trial(p,size,alpha,beta);
    calculate_probability();
    return function_probability;
}

double TrialFunction::get_probability_ratio(std::vector<Particle> p,int size,int move, double alpha, double beta){
    double val = 0;
    vec3 r;

    for(int i = 0; i<size;i++){
        if(i == move){
            r = p[i].next_r;
        }
        else{
            r = p[i].r;
        }

        val *= phi(r,alpha,beta);
    }


    return val/function_probability;
}



double TrialFunction::get_local_energy(int n, int dim){
    calculate_local_energy(n,dim);
    return local_energy;
}
