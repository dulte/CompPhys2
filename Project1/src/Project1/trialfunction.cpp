#include "trialfunction.h"
#include <iostream>

TrialFunction::TrialFunction(HarmonicOscillator *m_potential)
{
    potential = m_potential;
    beta=Parameters::beta;
    a=Parameters::a;
    dx=Parameters::dx;
    dimension = Parameters::dimension;
    N = Parameters::N;
    D = Parameters::D;

    h = 5e-8;

    allocate_empty_arrays();

    if(Parameters::a != 0){
        TrialFunction::f_func = &TrialFunction::f;
    }
    else{
        TrialFunction::f_func= &TrialFunction::f_id;
    }
    if(Parameters::D != 0){
        TrialFunction::greens_function_ratio_func = &TrialFunction::greens_function_ratio;
    }
    else{
        TrialFunction::greens_function_ratio_func= &TrialFunction::greens_function_ratio_id;
    }


}

void TrialFunction::allocate_empty_arrays(){

    std::vector<double> empty_N;
    for(int i=0; i<N;i++){
       empty_N.push_back(0.0);
    }
    for(int i=0; i<N;i++){
       distance_matrix.push_back(empty_N);
       distance_matrix_new.push_back(empty_N);
    }

    for(int i=0;i<N;i++){
        quantum_force_matrix.push_back(empty_N);
        quantum_force_matrix_new.push_back(empty_N);

    }
}

void TrialFunction::calculate_trial(std::vector<Particle> &p, double alpha)
{
    function_value = return_trial(p,alpha);
}

double TrialFunction::return_trial(std::vector<Particle> &p, double alpha)
{
    double val = 1;
    //std::vector<double> r;

    for(int i = 0; i<N;i++){
       // r = p[i].r;
        val *= phi(p[i].r,alpha)*(this->*f_func)(p);
    }

    return val;
}


void TrialFunction::calculate_probability(){
    function_probability = function_value*function_value;
    //function_probability_next_step = function_value_next_step*function_value_next_step;
}

void TrialFunction::calculate_local_energy(int n,int dim){
    local_energy = dim/2.0*n; // And external!!!
}

double TrialFunction::phi(std::vector<double> &r, double alpha)
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

double TrialFunction::f(std::vector<Particle> &p)
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

double TrialFunction::f_id(std::vector<Particle> &p){
    return 1.0;
}


void TrialFunction::quantum_force(std::vector<Particle> &p,double alpha)
{
    double dist=0;
    double dist_new=0;
    double grad_value=0;
    for(int i=0; i<N; i++){
        for(int j=0; j<i; j++){
            for(int k=0; k<dimension;k++){
                dist+=(p[i].r[k]-p[j].r[k])*(p[i].r[k]-p[j].r[k]);
            }
            dist=sqrt(dist);
            distance_matrix[i][j]=dist;
            distance_matrix[j][i]=dist;
            dist=0;
        }
    }

    for(int i=0; i<N;i++){
        for(int j=0; j<dimension;j++){
            for(int k=0; k<N;k++){
                if(k !=i){
                    if(distance_matrix[i][k]>a){
                        grad_value+=a*(p[i].r[j]-p[k].r[j])/(distance_matrix[i][k]*distance_matrix[i][k]*distance_matrix[i][k]);
                    }
                 }
            }
            if(j==2){
               quantum_force_matrix[i][j]=-4.0*alpha*beta*p[i].r[j]+grad_value;

            }
            else{
                quantum_force_matrix[i][j]=-4.0*alpha*p[i].r[j]+grad_value;
            }
            grad_value=0;
        }
    }
}

    void TrialFunction::quantum_force_new(std::vector<Particle> &p, double alpha, int chosen_particle)
    {
        double dist=0;
        double dist_new=0;
        double grad_value=0;
        std::vector<double> r;
        for(int i=0; i<N; i++){
            if(i == chosen_particle){
                r=p[i].next_r;
            }
            else{
                r=p[i].r;
            }
            for(int j=0; j<i; j++){
                for(int k=0; k<dimension;k++){
                    dist+=(r[k]-p[j].r[k])*(r[k]-p[j].r[k]);
                }
                dist=sqrt(dist);
                distance_matrix_new[i][j]=dist;
                distance_matrix_new[j][i]=dist;
                dist=0;
            }
        }

        for(int i=0; i<N;i++){
            if(i == chosen_particle){
                r=p[i].next_r;
            }
            else{
                r=p[i].r;
            }
            for(int j=0; j<dimension;j++){
                for(int k=0; k<N;k++){
                    if(k !=i){
                        if(distance_matrix[i][k]>a){
                            grad_value+=a*(r[j]-p[k].r[j])/(distance_matrix_new[i][k]*distance_matrix_new[i][k]*distance_matrix_new[i][k]);
                        }
                   }
                }
                if(j==2){
                   quantum_force_matrix_new[i][j]=-4.0*alpha*beta*r[j]+grad_value;

                }
                else{
                    quantum_force_matrix_new[i][j]=-4.0*alpha*r[j]+grad_value;
                }
                grad_value=0;
            }
        }
    }

double TrialFunction::greens_function_ratio(std::vector<Particle> &p, double alpha,int chosen_particle)
{
    double value=0;
    double value_new=0;
    double exponent = 0;
    double exponent_new = 0;
    double exponent_factor=D*dx;

    quantum_force(p, alpha);
    quantum_force_new(p, alpha,chosen_particle);
    std::vector<double> r;
    for(int i=0; i<N; i++){
        if(i==chosen_particle){
            r=p[i].next_r;
         }
         else{
             r=p[i].r;
         }
         for(int k=0;k<dimension;k++){
                exponent+=((p[i].r[k]-r[k])-exponent_factor*quantum_force_matrix[i][k])*((p[i].r[k]-r[k])-exponent_factor*quantum_force_matrix[i][k]);
                exponent_new+=((r[k]-p[i].r[k])-exponent_factor*quantum_force_matrix_new[i][k])*((r[k]-p[i].r[k])-exponent_factor*quantum_force_matrix_new[i][k]);

            }

        exponent=exponent/(4.0*exponent_factor);
        exponent_new=exponent_new/(4.0*exponent_factor);
        value+=exp(-exponent);
        value_new+=exp(-exponent_new);
        }
    return value_new/value;

}

double TrialFunction::greens_function_ratio_id(std::vector<Particle> &p, double alpha, int chosen_particle)
{
    return 1.0;
}



double TrialFunction::get_probability(std::vector<Particle> &p,int size,double alpha){
    calculate_trial(p,alpha);
    calculate_probability();
    return function_probability;
}


double TrialFunction::get_probability_ratio(std::vector<Particle> &p,int size,int move, double alpha){
    double val = 1;
    std::vector<double> r;
    double probability = get_probability(p,size,alpha);

    for(int i = 0; i<size;i++){
        if(i == move){
            r = p[i].next_r;
        }
        else{
            r = p[i].r;
        }

        val *= phi(r,alpha); //Can be optimized since this is a product of exp

    }
    return (this->*greens_function_ratio_func)(p, alpha, move)*val*val/probability;
}



double TrialFunction::get_local_energy(int n, int dim){
    calculate_local_energy(n,dim);
    return local_energy;
}

double TrialFunction::calculate_kinetic_energy(std::vector<Particle> &p,double alpha){
    double psi_minus = 0;
    double psi_plus = 0;

    double psi = return_trial(p,alpha);

    p_min = p;

    double kinetic_energy = 0;
    double potential_energy = 0;

    for(int i = 0; i<N;i++){
        potential_energy+=potential->get_external_potential(p[i].r);
        for(int j = 0; j<dimension;j++){
            p[i].r[j] += h;
            p_min[i].r[j]-= h;

            psi_minus = return_trial(p_min,alpha);
            psi_plus = return_trial(p,alpha);

            kinetic_energy -= (psi_plus+psi_minus - 2*psi);


            p[i].r[j] -= h;
            p_min[i].r[j]+= h;


        }

    }

    return (potential_energy+0.5*(kinetic_energy/(h*h))/psi);
}
