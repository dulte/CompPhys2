#include "simulation.h"
#include <iostream>
#include <math.h>

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

}

double Simulation::conjugate_gradient(double alpha_0, double b){
    double p_k=0;
    double B_k=b;
    double x_k=alpha_0;
    double x_k_1=0;
    double gradient_k=compute_local_energy_derivative(x_k);
    double local_energy_current=total_energy;
    double s_k=0;
    double gradient_k_1=0;
    double y_k=0;
    double alpha_k=0.1;
    double alpha_i=1;
    double r=1e-4;
    double temp=0;
    int i=0;
    int j=1;
    int max_iter=100;
    int max_iter_inner=10;
    int max_j=0;
    double tol=1e-8;

    double gradient_10_steps_ago = 1e6;



    while(i<max_iter){
        p_k=-gradient_k;//B_k;
        /*
        if(p_k < 0){
            p_k = -1;
        }
        else{
            p_k = 1;
        }*/
        std::cout << "Gradient: "<<gradient_k<<std::endl;
        /*
        while(j<max_iter_inner){
            temp=compute_local_energy_derivative(x_k+j*r);
            //std::cout<<total_energy << " "<<local_energy_current<<std::endl;
            if(total_energy<local_energy_current){
                alpha_k=j*r;
                std::cout<<"Found alpha k "<<alpha_k<<std::endl;
                break;
            }
            j++;
        }*/
        j=1;

        s_k=alpha_k*p_k;
        x_k_1=x_k+s_k;
        std::cout<< p_k << " " << s_k<<std::endl;
        //std::cout<<x_k_1<<std::endl;
        gradient_k_1=compute_local_energy_derivative(x_k_1);

        y_k=gradient_k_1-gradient_k;
        B_k=y_k/s_k;
        if(fabs(gradient_k) < tol || fabs(1-p_k/gradient_10_steps_ago) < tol){
            std::cout<<"MOM WE DID IT "<<gradient_k<<std::endl;
            break;
        }


        if(i%10 == 0){
            gradient_10_steps_ago = p_k;
        }
        gradient_k=gradient_k_1;
        local_energy_current=total_energy;
        x_k=x_k_1;
        i++;

    }



    return x_k;
}


void Simulation::check_derivative_of_energy(){

    DataDump<double> deriv_data("..//output//derivate_data.bin");
    DataDump<double> alpha_data("..//output//alpha_data.bin");

    double deriv = 0;

    for(double a = alpha_min; a < alpha_max; a+=alpha_step){
        deriv=compute_local_energy_derivative(a);
        deriv_data.push_back(deriv);
        alpha_data.push_back(a);
    }

    deriv_data.dump_all();
    alpha_data.dump_all();



}



/*
double Simulation::conjugate_gradient(double alpha_0, double b){
    double gradient=compute_local_energy_derivative(alpha_0);
    double alpha=alpha_0;
    double alpha_next=alpha_0-b*gradient;
    double alpha_prev=0;
    double gradient_prev=0;
    double tol=1e-14;
    int max_iter=300;
    int i=0;
    while(i<=max_iter){
        alpha=alpha_next;
        gradient=compute_local_energy_derivative(alpha);
        alpha_next=alpha-b*gradient;
        std::cout<<"Gradient "<<abs(gradient)<<std::endl;
        if(abs(gradient) < tol){
            std::cout<<"MOM WE DID IT"<<std::endl;
            break;
        }
        i++;
    }
    return alpha_next;

}
*/


double Simulation::compute_local_energy_derivative(double alpha){
    //std::cout<<alpha<<std::endl;

    energy = 0;
    total_energy = 0;
    energy_numerical=0;
    double local_energy_derivative=0;

    system->make_grid(alpha);
    for(int i = 0;i<MC_cycles;i++){
        energy = 0;
        for(int move = 0;move<N;move++){
            system->make_move_and_update(move);
            energy += system->check_acceptance_and_return_energy(move);
        }
            total_energy += energy/N;
        }
    total_energy =total_energy/MC_cycles;
    //std::cout<<"Derivative: " << system->expectation_derivative<<std::endl;
    //std::cout<<"D&E: " << system->expectation_derivative_energy<<std::endl;
    //std::cout<<"E: " << system->expectation_local_energy<<std::endl;

    double expectation_wavefunction_times_local_energy = system->expectation_derivative_energy/(N*MC_cycles);
    double expectation_wavefunction_times_expectation_local_energy = system->expectation_derivative/(N*MC_cycles)*total_energy;
    local_energy_derivative = 2.0*(expectation_wavefunction_times_local_energy-expectation_wavefunction_times_expectation_local_energy);
    //local_energy_derivative=(2.0/(N*MC_cycles))*(system->expectation_derivative_energy-(2.0/(N*MC_cycles))*system->expectation_derivative*system->expectation_local_energy);

    return local_energy_derivative;

}




void Simulation::run(int rank){

    //Initializes dump files

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

    //Main MC

    energy = 0;
    total_energy = 0;
    energy_numerical=0;

    for(double a = alpha_min; a < alpha_max; a+=alpha_step){
        system->make_grid(a);
        dump.push_back_stamp(a);
        std::cout << "a: " << a << std::endl;
        for(int i = 0;i<MC_cycles;i++){
            energy = 0;
            if(i%10000==0){
            //std::cout << "Mc steg: " << i << "/" << MC_cycles << std::endl;
            //std::cout << "Pos: " << *system->r << std::endl;
            }
            for(int move = 0;move<N;move++){
                system->make_move_and_update(move);
                energy += system->check_acceptance_and_return_energy(move);
            }
            //std::cout<<total_energy/N<<std::endl;
            total_energy += energy/N;
            dump.push_back(energy/N);
        }
        std::cout<<"Mean energy: "<<total_energy/(1.0*MC_cycles)<<std::endl;
        std::cout<<system->acceptance <<" "<<MC_cycles*N<<std::endl;
        energy_numerical=0;
        energy = 0;
        total_energy = 0;
        system->acceptance = 0;
    }
    dump.dump_all();
}


