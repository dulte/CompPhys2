#include "simulation.h"
#include <iostream>

Simulation::Simulation(System *m_system)
{
    system = m_system;
}

double Simulation::conjugate_gradient(double alpha_0, double b){
    double p_k=0;
    double B_k=b;
    double x_k=alpha_0;
    double x_k_prev=alpha_0;
    double gradient_k=compute_local_energy_derivative(x_k);
    double local_energy_current=total_energy;
    double s_k=0;
    double gradient_k_prev=0;
    double y_k=0;
    double alpha_k=0.01;
    double alpha_i=1;
    double r=1e-4;
    double temp=0;
    int i=0;
    int j=1;
    int max_iter=500;
    int max_iter_inner=10;
    int max_j=0;
    double tol=1e-5;

    double gradient_10_steps_ago = 1e6;

    gradient_k_prev=compute_local_energy_derivative(x_k);
    //x_k += alpha_k*gradient_k_prev;
    p_k = -gradient_k_prev;


    while(i<max_iter){


        std::cout << "Gradient: "<<gradient_k<<std::endl;
        s_k=alpha_k*p_k;
        x_k = x_k_prev + s_k;

        if(fabs(x_k)>1e5 || fabs(gradient_k)>1e10 || gradient_k != gradient_k){
            std::cout << "Reset" << std::endl;
            x_k = alpha_0;
            gradient_k_prev=compute_local_energy_derivative(x_k);
            p_k = -gradient_k_prev;
            alpha_k = 0.01;
            x_k += p_k*alpha_k;
            gradient_k=compute_local_energy_derivative(x_k);
            local_energy_current = 0;
            i = 0;
        }



        gradient_k=compute_local_energy_derivative(x_k);


        if(fabs(gradient_k - gradient_k_prev) > 1e-10){
            alpha_k = (x_k - x_k_prev)/(gradient_k - gradient_k_prev);
        }
        else{
            std::cout << "Eeeeh" << std::endl;
            alpha_k = 0.001;
        }

        p_k = -gradient_k;

        //gradient_k_1=compute_local_energy_derivative(x_k_1);

        //y_k=gradient_k_1-gradient_k;
        //B_k=y_k/s_k;


        if(fabs(gradient_k) < tol || fabs(1-p_k/gradient_10_steps_ago) < tol){
            std::cout<<"MOM WE DID IT "<<gradient_k<<std::endl;
            break;
        }


        if(i%10 == 0){
            gradient_10_steps_ago = p_k;
        }
        gradient_k_prev=gradient_k;
        local_energy_current=total_energy;
        x_k_prev=x_k;
        i++;

    }



    return x_k;
}

double Simulation::compute_local_energy_derivative(double alpha){
    energy = 0;
    total_energy = 0;
    int move = 0;
    double local_energy_derivative=0;
    int fast_MC_cycles = static_cast<int>(MC_cycles/100.);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    system->make_grid(alpha);
    for(int i = 0;i<fast_MC_cycles;i++){
        energy = 0;
        move = (int)distribution(gen);
        system->make_move_and_update(move);
        energy += system->check_acceptance_and_return_energy(move);

        total_energy += energy;
    }

    total_energy =total_energy/fast_MC_cycles;
    double expectation_wavefunction_times_local_energy = system->expectation_derivative_energy/(N*fast_MC_cycles);
    double expectation_wavefunction_times_expectation_local_energy = system->expectation_derivative/(N*fast_MC_cycles)*total_energy;
    local_energy_derivative = 2.0*(expectation_wavefunction_times_local_energy-expectation_wavefunction_times_expectation_local_energy);
    //local_energy_derivative=(2.0/(N*MC_cycles))*(system->expectation_derivative_energy-(2.0/(N*MC_cycles))*system->expectation_derivative*system->expectation_local_energy);

    return local_energy_derivative;

}

void Simulation::data_for_derivated(){
    DataDump<double> deriv_data("..//output//deriv_data.bin");
    DataDump<double> data("..//output//data.bin");
    DataDump<double> alphas("..//output//alphas.bin");
    for(double a = alpha_min; a < alpha_max; a+=alpha_step){
        std::cout << "a: " << a << std::endl;
        alphas.push_back(a);
        double deriv = compute_local_energy_derivative(a);
        data.push_back(total_energy);
        deriv_data.push_back(deriv);
    }
    data.dump_all();
    alphas.dump_all();
    deriv_data.dump_all();

}




void Simulation::initiate(){
    alpha_max = Parameters::alpha_max;
    alpha_min = Parameters::alpha_min;
    alpha_step = (alpha_max - alpha_min)/(double)Parameters::alpha_num;
    MC_cycles = Parameters::MC_cycles;
    N = Parameters::N;
}


void Simulation::run(int rank){

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    energy = 0;
    double total_energy;
    int move = 0;

    std::string filename = "..//output//data";
    filename.append(std::to_string(rank));
    filename.append(".bin");

    std::string stampname = "..//output//stamp";
    stampname.append(std::to_string(rank));
    stampname.append(".bin");

    DataDump<double> dump(filename,stampname);

    DataDump<double> data("..//output//data.bin");
    if(rank == 0){
         dump.dump_metadata("..//output//metadata.txt");

    }



    for(double a = alpha_min; a < alpha_max; a+=alpha_step){

        system->make_grid(a);
        dump.push_back_stamp(a);
        std::cout << "a: " << a << std::endl;
        for(int i = 0;i<MC_cycles;i++){
            energy = 0;
            move = (int)distribution(gen);
            system->make_move_and_update(move);
            energy += system->check_acceptance_and_return_energy(move);
            dump.push_back(energy);
            total_energy += energy;

        }
        //dump.push_back(energy/(MC_cycles));
        if(rank == 0){
            data.push_back(total_energy/(MC_cycles));
        }
        std::cout << "Energy " << total_energy/(MC_cycles) << std::endl;
        total_energy = 0;
    }
    dump.dump_all();
    if(rank == 0){
        data.dump_all();
    }

}


void Simulation::run(int rank,double alpha){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    energy = 0;
    double total_energy;
    int move = 0;

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
        move = (int)distribution(gen);
        system->make_move_and_update(move);
        energy += system->check_acceptance_and_return_energy(move);
        dump.push_back(energy);
        total_energy += energy;

    }
    std::cout << "Energy " << total_energy/(MC_cycles) << std::endl;
    dump.dump_all();
}


void Simulation::oneBodyDensity(double optimal_alpha, double r_step,double r_min, double r_max){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0,N);

    int move = 1;
    double pi = 3.14159265359;




    DataDump<double> r_packet("..//output//r_positions.bin");

    DataDump<std::vector<double>> density_packet("..//output//density.bin");





    int r_num = (int)(r_max - r_min)/(r_step);

    std::cout << r_num << " " << (r_max - r_min)/(r_step) << std::endl;

    std::vector<double> density;
    std::vector<double> rs;
    std::vector<double> volume;
    for(double r = r_min;r<=r_max;r+=r_step){
        rs.push_back(r);
        density.push_back(0);

       // for(int i = 0; i<r_num;i++){


        if(Parameters::dimension ==1){
            volume.push_back(r+r_step-r);
        }
        else if(Parameters::dimension == 2){
            volume.push_back(M_PI*((r+r_step)*(r+r_step) - r*r));
        }
        else{
            volume.push_back(4./3*M_PI*((r+r_step)*(r+r_step)*(r+r_step) - r*r*r));
        }
    //}

    }






    r_packet.push_back(r_min);
    r_packet.push_back(r_max);
    r_packet.push_back((double)r_num);
    r_packet.push_back(r_step);


    double distance_from_origo = 0;
    double wave_function_squared = 0;
    int bin = 0;


    system->make_grid(optimal_alpha);

    for(int i = 0;i<MC_cycles;i++){
        move = (int)distribution(gen);
        system->make_move_and_update(move);
        system->check_acceptance_and_return_energy(move);
        distance_from_origo = system->r.col(move).norm();

        if(distance_from_origo > r_max){
            std::cout << "Calculated distance is larger that r_max" << std::endl;
            exit(1);
        }

        bin = (int)r_num*distance_from_origo/r_max;

        density[bin] += 1./(MC_cycles*volume[bin]);
    }

    density_packet.push_back(density);


    r_packet.dump_all();
    density_packet.dump_all();

}
