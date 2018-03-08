#include "simulation.h"
#include <iostream>

Simulation::Simulation(System *m_system)
{
    system = m_system;
}

void Simulation::initiate(){
    alpha_max = Parameters::alpha_max;
    alpha_min = Parameters::alpha_min;
    alpha_step = (alpha_max - alpha_min)/(double)Parameters::alpha_num;
    MC_cycles = Parameters::MC_cycles;
    N = Parameters::N;

    //system->make_grid(alpha_min);
}

double Simulation::conjugate_gradient(double alpha_0, double b){
    double gradient=compute_local_energy_derivative(alpha_0);
    double alpha=alpha_0;
    double alpha_next=alpha_0-b*gradient;
    double alpha_prev=0;
    double gradient_prev=0;
    double tol=1e-14;
    int max_iter=30;
    int i=0;
    while(i<=max_iter){
        alpha=alpha_next;
        gradient=compute_local_energy_derivative(alpha);
        std::cout<<"Gradient "<<gradient<<std::endl;
        alpha_next=alpha-b*gradient;
        if(abs(gradient) < tol){
            std::cout<<"MOM WE DID IT"<<std::endl;
            break;
        }
        i++;
    }
    return alpha_next;

}


double Simulation::compute_local_energy_derivative(double alpha){
    energy = 0;

    system->make_grid(alpha);
    std::cout << "a: " << alpha << std::endl;
    for(int i = 0;i<MC_cycles;i++){
        for(int move = 0;move<N;move++){
            system->make_move_and_update(move);
            energy += system->check_acceptance_and_return_energy(move);
        }
    }
    //std::cout<<"Derivative: " << system->expectation_derivative<<std::endl;
    //std::cout<<"D&E: " << system->expectation_derivative_energy<<std::endl;
    //std::cout<<"E: " << system->expectation_local_energy<<std::endl;
    return (2.0/(N*MC_cycles))*(system->expectation_derivative_energy-system->expectation_derivative*system->expectation_local_energy);


}


void Simulation::run(){

    energy = 0;
    energy_numerical=0;

    for(double a = alpha_min; a < alpha_max; a+=alpha_step){
        system->make_grid(a);
        std::cout << "a: " << a << std::endl;
        for(int i = 0;i<MC_cycles;i++){
            if(i%10000==0){
            //std::cout << "Mc steg: " << i << "/" << MC_cycles << std::endl;
            //std::cout << "Pos: " << *system->r << std::endl;
            }
            for(int move = 0;move<N;move++){
                system->make_move_and_update(move);
                energy += system->check_acceptance_and_return_energy(move);
            }
            energy_numerical+=system->calculate_energy();
        }
    std::cout<<"Mean energy: "<<energy_numerical/(1.0*MC_cycles)<<std::endl;
    std::cout<<system->acceptance <<" "<<MC_cycles*N<<std::endl;
    energy_numerical=0;
    system->acceptance = 0;
    }
}


