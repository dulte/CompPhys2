#include "simulation.h"
#include <iostream>


Simulation::Simulation(std::shared_ptr<System> m_system)
{
    m_system->update_alpha(0);
    system = m_system;
}

void Simulation::initiate(int m_size, int m_alpha_min,int m_alpha_max,int m_alpha_num){
    size = m_size;
    alpha_max = m_alpha_max;
    alpha_min = m_alpha_min;
    alpha_step = (m_alpha_max - m_alpha_min)/(double)m_alpha_num;

    system->grid_setup(size,alpha_min);
}

void Simulation::run(int m_MCsteps){

    double energy = 0;
    double energy_squared = 0;
    double delta_error = 0;


    for(double a = alpha_min; a <= alpha_max; a = a + alpha_step){

        for (int i = 0;i<m_MCsteps;i++){

            system->update_alpha(a);
            system->propose_step();

            delta_error = system->check_acceptance_and_return_energy();

            energy = energy + delta_error;
            energy_squared = energy_squared + delta_error*delta_error;

        }
    }

}



