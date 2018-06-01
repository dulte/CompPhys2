#include "simulation.h"
#include <iostream>

Simulation::Simulation(System *m_system)
{
    system = m_system;
    MC_cycles = Parameters::MC_cycles;

}



Eigen::ArrayXd Simulation::stochastic_descent(Eigen::ArrayXd x_0){
    int max_iter = 100;
    int i = 0;
    double A = 20;
    Eigen::ArrayXd t = Eigen::ArrayXd::Ones(x_0.size())*A;
    Eigen::ArrayXd x = x_0;
    Eigen::ArrayXd x_prev = x_0;
    Eigen::ArrayXd gradient = Eigen::ArrayXd::Zero(x_0.size());
    double tol = 1e-5;

    std::string gradient_filename = "..//output//gradient_data_";
    gradient_filename.append(std::to_string(Parameters::N));
    gradient_filename.append("_");
    gradient_filename.append(std::to_string(Parameters::learning_rate));

    std::string energy_filename = "..//output//energy_data_";
    energy_filename.append(std::to_string(Parameters::N));
    energy_filename.append("_");
    energy_filename.append(std::to_string(Parameters::learning_rate));

    DataDump<double> gradient_dump(gradient_filename);
    DataDump<double> energy_dump(energy_filename);





    while(i < max_iter){
        calculate_gradient(x,gradient);
        //std::cout << "Gradient: " << (gradient) << std::endl;
        x = x_prev - Parameters::learning_rate*gradient;
        x_prev = x;

        gradient_dump.push_back(((Eigen::VectorXd)gradient).squaredNorm());
        energy_dump.push_back(total_energy);



        if(((Eigen::VectorXd)gradient).squaredNorm() < tol){
            std::cout << "Mom, we did it" << std::endl;
            std::cout << "Norm of gradient: " << ((Eigen::VectorXd)gradient).squaredNorm() << std::endl;
            std::cout << "Energy: " << total_energy << std::endl;
            std::cout << "Number of iterations: " << i+1 << std::endl;
            break;
        }
        //std::cout << gradient << std::endl;
        std::cout << "Norm of gradient: " << ((Eigen::VectorXd)gradient).squaredNorm() << std::endl;
        i++;

    }

    //Uncomment this to write energy and gradient data to file
    //gradient_dump.dump_all();
    //energy_dump.dump_all();
    return x;
}


void Simulation::calculate_gradient(Eigen::ArrayXd &x,Eigen::ArrayXd &gradient){

    int move = 0;
    double local_energy_derivative=0;
    double local_energy = 0;
    total_energy = 0;

    double variable_derivative = 0;
    int total_size = x.size();
    int M = Parameters::P*Parameters::dimension;



    Eigen::ArrayXd E_L_times_derivatives = Eigen::ArrayXd::Zero(total_size);
    Eigen::ArrayXd derivatives = Eigen::ArrayXd::Zero(total_size);

    //This lowers the number of MC step by a factor 100.
    //This is done because we dont need as many steps, and to make this go faster.
    int fast_MC_cycles = static_cast<int>(MC_cycles/10.);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,Parameters::P);
    //A simple simulation for the given alpha
    system->make_grid(x);

    for(int i = 0;i<fast_MC_cycles;i++){

        if(Parameters::gibbs){

            local_energy = system->gibbs_sample_and_return_energy();
        }else{
            move = static_cast<int>(distribution(gen));
            system->make_move_and_update(move);
            local_energy = system->check_acceptance_and_return_energy(move);
        }

        total_energy += local_energy;

        for(int k = 0;k<total_size;k++){
            if(k<M){
                variable_derivative = system->d_psi_da(k);

            }else if(k>=M && k<(Parameters::N+M)){
                variable_derivative = system->d_psi_db(k-M);
            }else{
                int w_index = k-(M+Parameters::N);
                int column = w_index/M;
                int row = w_index%M;
                variable_derivative = system->d_psi_dw(row,column);

            }

            E_L_times_derivatives(k) += variable_derivative*local_energy;
            derivatives(k) += variable_derivative;

            if(variable_derivative != variable_derivative){
                exit(1);
            }


        }

    }

    E_L_times_derivatives /= fast_MC_cycles;
    derivatives /= fast_MC_cycles;
    total_energy /= fast_MC_cycles;

    std::cout << "total energy: " << total_energy << std::endl;
    std::cout << system->number_accept/((double)fast_MC_cycles) << std::endl;
    /*std::cout << "#######" << std::endl;
    std::cout << total_energy*derivatives << std::endl;*/
    std::cout << "-----------------" << std::endl;

    gradient = 2*(E_L_times_derivatives - total_energy*derivatives);
    std::cout << "Gradients: " << std::endl;
    for(int k = 0;k<total_size;k++){
        /*if(k<M){
            std::cout << "a: " << gradient[k] << std::endl;

        }else if(k>=M && k<(Parameters::N+M)){
            std::cout << "b: " << gradient[k-M] << std::endl;
        }else{

            std::cout << "w: " << gradient[k - (M+N)] << std::endl;
        }

        E_L_times_derivatives(k) += variable_derivative*local_energy;
        derivatives(k) += variable_derivative;

        if(variable_derivative != variable_derivative){
            exit(1);
        }

        */
    }

    std::cout << "_________________" << std::endl;


}




//Does the main Metropolise Simulation over the alpha interval given in the parameter file.
void Simulation::run(int rank,Eigen::ArrayXd &x){

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,Parameters::P);

    energy = 0;
    double total_energy;
    int move = 0;

    //Makes filenames containing the rank of the process
    std::string filename = "..//output//data";
    //filename.append(std::to_string(rank));
    //filename.append(std::to_string(Parameters::dx));
    filename.append(".bin");


    //First filename holdes the main data(energy), while the second filename holdes the stamp(Alphas)
    DataDump<double> dump(filename);


    std::string accept_filename = "..//output//accept_data";
    //filename.append(std::to_string(rank));
    //accept_filename.append(std::to_string(Parameters::dx));
    accept_filename.append(".bin");

    DataDump<double> accept_dump(accept_filename);


    if(rank == 0){
         dump.dump_metadata("..//output//metadata.txt");

    }


    //Makes a new system with the given alpha
    system->make_grid(x);

    std::cout<<MC_cycles << std::endl;
    int count = 0;

    //The main MC loop
    for(int i = 0;i<MC_cycles;i++){
        energy = 0;
        if (!Parameters::gibbs){
            move = static_cast<int>(distribution(gen));

            //Makes the move and updates the relevant data
            system->make_move_and_update(move);

            //Checks if the move is accepted, and returns the local energy.
            energy += system->check_acceptance_and_return_energy(move);
        }else{
            energy += system->gibbs_sample_and_return_energy();
        }
        dump.push_back(energy);
        accept_dump.push_back(system->number_accept/double(i+1));
        total_energy += energy;

        count++;

    }

    std::cout << "When done: " << count << std::endl;

    std::cout << "Energy " << total_energy/(MC_cycles) << std::endl;
    std::cout << "Acceptance rate: " << (double)system->number_accept/double(MC_cycles) << std::endl;
    total_energy = 0;

    dump.dump_all();
    accept_dump.dump_all();

}


