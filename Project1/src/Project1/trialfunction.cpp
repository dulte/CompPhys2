#include "trialfunction.h"
#include <iostream>

TrialFunction::TrialFunction(std::shared_ptr<Potential> m_potential)
{
    potential = m_potential;
    beta=Parameters::beta;
    a=Parameters::a;
    dx=Parameters::dx;
    if(Parameters::a != 0){
        TrialFunction::f_func = &TrialFunction::f;
    }
    else{
        TrialFunction::f_func= &TrialFunction::f_id;
    }

}

void TrialFunction::calculate_trial(std::vector<Particle> p, int size, double alpha)
{
    double val = 1;
    std::vector<double> r;

    for(int i = 0; i<size;i++){
        r = p[i].r;
        val *= phi(r,alpha)*(this->*f_func)(p);
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

double TrialFunction::phi(std::vector<double> r, double alpha)
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

double TrialFunction::f(std::vector<Particle> p)
{
    int size = p.size();
    double dist=0;
    double val = 0;
    for(int i = 0;i<size;i++){
        for(int j = 0; j<i;j++){
            for(int k = 0; k<3; k++){
                dist+=(p[i].r[k]-p[j].r[k])*(p[i].r[k]-p[j].r[k]);
            }
            dist=sqrt(dist);
            if(dist>a){
                val*=1.0-a/dist;
            }
       }
    }

    return val;
}

double TrialFunction::f_id(std::vector<Particle> p){
    return 1.0;
}



double TrialFunction::get_probability(std::vector<Particle> p,int size,double alpha){
    calculate_trial(p,size,alpha);
    calculate_probability();
    return function_probability;
}


double TrialFunction::get_probability_ratio(std::vector<Particle> p,int size,int move, double alpha){
    double val = 1;
    std::vector<double> r;

    for(int i = 0; i<size;i++){
        if(i == move){
            r = p[i].next_r;
        }
        else{
            r = p[i].r;
        }

        val *= phi(r,alpha); //Can be optimized since this is a product of exp

    }
    return val*val/function_probability;
}



double TrialFunction::get_local_energy(int n, int dim){
    calculate_local_energy(n,dim);
    return local_energy;
}
