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


void Simulation::run(){

    energy = 0;
    energy_numerical=0;

    for(double a = alpha_min; a < alpha_max; a+=alpha_step){
        system->make_grid(a);
        std::cout << "a: " << a << std::endl;
        for(int i = 0;i<MC_cycles;i++){
            if(i%100==0){
            std::cout << "Mc steg: " << i << "/" << MC_cycles << std::endl;
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
