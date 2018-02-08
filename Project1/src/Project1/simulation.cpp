#include "simulation.h"
#include <iostream>


Simulation::Simulation(std::shared_ptr<System> m_system)
{
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

    DataDump<double> dump("..//output//data.txt");
    DataDump<std::vector<double>> position_dump("..//output//positiondata.txt");

    double energy = 0;
    double energy_squared = 0;
    double delta_energy = 0;
    double kinetic_energy = 0;


    for(double a = alpha_min; a <= alpha_max; a = a + alpha_step){

        for (int i = 0;i<m_MCsteps;i++){

            system->update_alpha(a);
            system->propose_step();

            delta_energy = system->check_acceptance_and_return_energy();
            kinetic_energy = system->trial_function->calculate_kinetic_energy(system->particles,a);

            energy = energy + delta_energy;
            energy_squared = energy_squared + delta_energy*delta_energy;
            std::cout << "Energy " << kinetic_energy << std::endl;
            dump.push_back(energy);
            position_dump.push_back(system->get_postions());

        }
    }
    dump.dump_all();
    position_dump.dump_all();

}



