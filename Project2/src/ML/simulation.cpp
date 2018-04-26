#include "simulation.h"
#include <iostream>

Simulation::Simulation(System *m_system)
{
    system = m_system;
    MC_cycles = Parameters::MC_cycles;

}

inline Eigen::ArrayXd f(const Eigen::ArrayXd & X){
    double f_max = 1;
    double f_min = -1;
    double f_ratio = f_max/f_min;
    double omega = 1;

    return  (f_min - (f_max-f_min)*Eigen::inverse(1-f_ratio*Eigen::exp(-X/omega)));
}


inline Eigen::ArrayXd step_length(Eigen::ArrayXd & X, const int & A,Eigen::ArrayXd & t){
    int a = 1;
    Eigen::ArrayXd f_val = f(X);
    Eigen::ArrayXd delta = a*Eigen::inverse(t+A);
    t = (t+f_val).cwiseMax(0);
    return delta;
}


Eigen::ArrayXd Simulation::stochastic_descent(Eigen::ArrayXd x_0){
    int max_iter = 2000;
    int i = 0;
    double A = 2;
    Eigen::ArrayXd  t = Eigen::ArrayXd::Ones(x_0.size())*A;
    Eigen::ArrayXd x = x_0;
    Eigen::ArrayXd x_prev = x_0;
    Eigen::ArrayXd gradient = Eigen::ArrayXd::Zero(x_0.size());



    while(i < max_iter){
        calculate_gradient(x,gradient);
        //std::cout << "Gradient: " << gradient << std::endl;
        x = x_prev - 0.01*gradient;//step_length(x_prev,A,t)*gradient;
        x_prev = x;
        i++;

    }
    return x;
}


void Simulation::calculate_gradient(Eigen::ArrayXd &x,Eigen::ArrayXd &gradient){

    int move = 0;
    double local_energy_derivative=0;
    double local_energy = 0;
    double total_energy = 0;

    double variable_derivative = 0;
    int total_size = x.size();
    int M = Parameters::P*Parameters::dimension;

    Eigen::ArrayXd E_L_times_derivatives = Eigen::ArrayXd::Zero(total_size);
    Eigen::ArrayXd derivatives = Eigen::ArrayXd::Zero(total_size);

    //This lowers the number of MC step by a factor 100.
    //This is done because we dont need as many steps, and to make this go faster.
    int fast_MC_cycles = static_cast<int>(MC_cycles/100.);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    //A simple simulation for the given alpha
    system->make_grid(x);

    for(int i = 0;i<fast_MC_cycles;i++){
        move = static_cast<int>(distribution(gen));
        system->make_move_and_update(move);

        local_energy = system->check_acceptance_and_return_energy(move);
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
    /*std::cout << "#######" << std::endl;
    std::cout << total_energy*derivatives << std::endl;*/
    std::cout << "-----------------" << std::endl;

    gradient = 2*(E_L_times_derivatives - total_energy*derivatives);

}







//Runs a simulation for a given alpha, and returns the derivative of the local energy, needed for gradien descent
double Simulation::compute_local_energy_derivative(double alpha){
    energy = 0;
    total_energy = 0;
    int move = 0;
    double local_energy_derivative=0;

    //This lowers the number of MC step by a factor 100.
    //This is done because we dont need as many steps, and to make this go faster.
    int fast_MC_cycles = static_cast<int>(MC_cycles/100.);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    //A simple simulation for the given alpha
    system->make_grid(alpha);
    for(int i = 0;i<fast_MC_cycles;i++){
        energy = 0;
        move = static_cast<int>(distribution(gen));
        system->make_move_and_update(move);
        energy += system->check_acceptance_and_return_energy(move);

        total_energy += energy;
    }

    total_energy =total_energy/fast_MC_cycles;
    double expectation_wavefunction_times_local_energy = system->expectation_derivative_energy/(N*fast_MC_cycles);
    double expectation_wavefunction_times_expectation_local_energy = system->expectation_derivative/(N*fast_MC_cycles)*total_energy;
    local_energy_derivative = 2.0*(expectation_wavefunction_times_local_energy-expectation_wavefunction_times_expectation_local_energy);


    return local_energy_derivative;

}



void Simulation::initiate(){

}




//Does the main Metropolise Simulation over the alpha interval given in the parameter file.
void Simulation::run(int rank){

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    energy = 0;
    double total_energy;
    int move = 0;

    //Makes filenames containing the rank of the process
    std::string filename = "..//output//data";
    filename.append(std::to_string(rank));
    filename.append(".bin");

    std::string stampname = "..//output//stamp";
    stampname.append(std::to_string(rank));
    stampname.append(".bin");

    //First filename holdes the main data(energy), while the second filename holdes the stamp(Alphas)
    DataDump<double> dump(filename,stampname);


    if(rank == 0){
         dump.dump_metadata("..//output//metadata.txt");

    }



    for(double a = alpha_min; a < alpha_max; a+=alpha_step){


        //Makes a new system with the given alpha
        system->make_grid(a);
        dump.push_back_stamp(a);
        std::cout << "a: " << a << std::endl;


        //The main MC loop
        for(int i = 0;i<MC_cycles;i++){
            energy = 0;
            move = static_cast<int>(distribution(gen));

            //Makes the move and updates the relevant data
            system->make_move_and_update(move);

            //Checks if the move is accepted, and returns the local energy.
            energy += system->check_acceptance_and_return_energy(move);

            dump.push_back(energy);
            total_energy += energy;

        }

        std::cout << "Energy " << total_energy/(MC_cycles) << std::endl;
        std::cout << "Acceptance rate: " << (double)system->number_accept/double(MC_cycles) << std::endl;
        total_energy = 0;
    }
    dump.dump_all();

}



//Same as the run above, but does the simulation for a given alpha
void Simulation::run(int rank,double alpha){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    energy = 0;
    double total_energy;
    int move = 0;

    //Makes filenames containing the rank of the process
    std::string filename = "..//output//data";
    filename.append(std::to_string(rank));
    filename.append(".bin");

    std::string stampname = "..//output//stamp";
    stampname.append(std::to_string(rank));
    stampname.append(".bin");

    DataDump<double> dump(filename,stampname);

    if(rank == 0){
         dump.dump_metadata("..//output//metadata.txt");

    }


    system->make_grid(alpha);
    dump.push_back_stamp(alpha);

    for(int i = 0;i<MC_cycles;i++){
        energy = 0;
        move = static_cast<int>(distribution(gen));
        system->make_move_and_update(move);
        energy += system->check_acceptance_and_return_energy(move);
        dump.push_back(energy);
        total_energy += energy;

    }
    std::cout << "Energy " << total_energy/(MC_cycles) << std::endl;
    dump.dump_all();
}

/*
 This function calculates the onebody density by simulating the system, and
 seeing at which raduis from the origin the newly moved particled is. In the python script,
 the number of particles is ajusted for.
*/
void Simulation::oneBodyDensity(double optimal_alpha, double r_step,double r_min, double r_max){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    int move = 1;
    double distance_from_origo = 0;
    int bin = 0;
    int r_num = (int)(r_max - r_min)/(r_step);

    std::vector<double> density;
    std::vector<double> rs;
    std::vector<double> volume;

    DataDump<double> r_packet("..//output//r_positions.bin");
    DataDump<std::vector<double>> density_packet("..//output//density.bin");
    DataDump<double> volume_factor_packet("..//output//volume.bin");




    for(double r = r_min;r<=r_max;r+=r_step){
        rs.push_back(r);
        density.push_back(0);


        //Makes the volume element
        volume.push_back(4*M_PI*(r+r_step)*(r+r_step)*r_step);
        volume_factor_packet.push_back(4*M_PI*(r+r_step)*(r+r_step)*r_step);
    }


    r_packet.push_back(r_min);
    r_packet.push_back(r_max);
    r_packet.push_back((double)r_num);
    r_packet.push_back(r_step);


    system->make_grid(optimal_alpha);

    for(int i = 0;i<MC_cycles;i++){

        //Makes a move, and check its acceptance
        move = static_cast<int>(distribution(gen));
        system->make_move_and_update(move);
        system->check_acceptance_and_return_energy(move);

        //Calculates the distance from the origo
        distance_from_origo = system->r.col(move).norm();

        if(distance_from_origo > r_max){
            std::cout << "Calculated distance is larger that r_max" << std::endl;
            exit(1);
        }

        //Calulates the bin where the particles is, and updates the density array
        bin = (int)r_num*distance_from_origo/r_max;
        density[bin] += 1./(MC_cycles*volume[bin]);

    }

    density_packet.push_back(density);


    r_packet.dump_all();
    density_packet.dump_all();
    volume_factor_packet.dump_all();
}
