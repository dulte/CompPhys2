#include "simulation.h"
#include <iostream>


Simulation::Simulation(System *m_system)
{
    system = m_system;
}

void Simulation::initiate(){
    size = Parameters::N;
    alpha_max = Parameters::alpha_max;
    alpha_min = Parameters::alpha_min;
    alpha_step = (alpha_max - alpha_min)/(double)Parameters::alpha_num;

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
        system->grid_setup(size,a);
        system->update_alpha(a);
        for (int i = 0;i<m_MCsteps;i++){
            system->propose_step();
            delta_energy = system->check_acceptance_and_return_energy();
            //kinetic_energy = system->trial_function->calculate_kinetic_energy(system->particles,a);
            energy += kinetic_energy;
        }

        std::cout << "Energy " << energy/(double)m_MCsteps << std::endl;
        energy = 0;
    }

    dump.dump_all();
    position_dump.dump_all();

}



