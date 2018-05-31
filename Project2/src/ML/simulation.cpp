#include "simulation.h"
#include <iostream>

Simulation::Simulation(System *m_system)
{
    system = m_system;
    MC_cycles = Parameters::MC_cycles;

}

inline Eigen::ArrayXd f(const Eigen::ArrayXd & X){
    double f_max = 1.;
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
    int max_iter = 100;
    int i = 0;
    double A = 20;
    Eigen::ArrayXd t = Eigen::ArrayXd::Ones(x_0.size())*A;
    Eigen::ArrayXd x = x_0;
    Eigen::ArrayXd x_prev = x_0;
    Eigen::ArrayXd gradient = Eigen::ArrayXd::Zero(x_0.size());
    double tol = 1e-3;

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
        x = x_prev - Parameters::learning_rate*gradient;//step_length(x_prev,A,t)*gradient;
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

    gradient_dump.dump_all();
    energy_dump.dump_all();
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







//Runs a simulation for a given alpha, and returns the derivative of the local energy, needed for gradien descent
double Simulation::compute_local_energy_derivative(double alpha){
    energy = 0;
    total_energy = 0;
    int move = 0;
    double local_energy_derivative=0;

    //This lowers the number of MC step by a factor 100.
    //This is done because we dont need as many steps, and to make this go faster.
    int fast_MC_cycles = static_cast<int>(MC_cycles/10.);

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
void Simulation::run(int rank,Eigen::ArrayXd &x){

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,Parameters::P);

    energy = 0;
    double total_energy;
    int move = 0;

    //Makes filenames containing the rank of the process
    std::string filename = "..//output//data_";
    //filename.append(std::to_string(rank));
    filename.append(std::to_string(Parameters::dx));
    filename.append(".bin");


    //First filename holdes the main data(energy), while the second filename holdes the stamp(Alphas)
    DataDump<double> dump(filename);


    std::string accept_filename = "..//output//accept_data_";
    //filename.append(std::to_string(rank));
    accept_filename.append(std::to_string(Parameters::dx));
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





/*
 This function calculates the onebody density by simulating the system, and
 seeing at which raduis from the origin the newly moved particled is. In the python script,
 the number of particles is ajusted for.
*/
void Simulation::oneBodyDensity(double optimal_alpha, double r_step,double r_min, double r_max){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,Parameters::P);

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
