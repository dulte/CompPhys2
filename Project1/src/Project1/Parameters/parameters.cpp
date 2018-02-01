#include "parameters.h"


void Parameters::read_parameters(std::string location){
    std::ifstream infile;
    std::string line;

    infile.open(location);
    if(infile.fail()){
        std::cout << "Could not find parameter file at " << location << std::endl;
        exit(EXIT_FAILURE);
    }
    infile.clear();
    infile.seekg(0,std::ios::beg);
    for (std::string l; getline(infile,l);)
    {
        std::stringstream ss(l);

        std::string name;
        double variable;

        ss >> name >> variable;
        //std::cout << "Ting: " << name[0] << std::endl;

        if (name.front() == '#'){
            continue;
        }
        else if(name == "text"){
            text = variable;
        }
        else if(name == "MC_cycles"){
            MC_cycles = variable;
            if(MC_cycles>0){
                MC_cycles_set=true;
            }
        }
        else if(name == "alpha_min"){
            alpha_min = variable;
            alpha_min_set=true;
        }
        else if(name == "alpha_max"){
            alpha_max = variable;
            alpha_max_set=true;
        }
        else if(name == "alpha_num"){
            alpha_num = variable;
            alpha_num_set=true;
        }

        else if(name == "beta"){
            beta = variable;
            beta_set=true;
        }
        else if(name == "omega"){
            omega = variable;
            omega_set=true;
        }

        else if(name == "omega_z"){
            omega_z=variable;
            omega_z_set=true;
        }

        else if(name == "a"){
            a=variable;
            a_set=true;
        }

        else{
            std::cout << "Unknonw Variable found: " << name << std::endl;
            exit(EXIT_FAILURE);
        }
    }


    if(!alpha_max_set){
        std::cout << "Max alpha not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!alpha_min_set){
        std::cout << "Min alpha not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!alpha_num_set){
        std::cout << "Numbers of alpha not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!beta_set){
        std::cout << "Beta not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!omega_set){
        std::cout << "Omega not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!omega_z_set){
        std::cout << "Omega_z not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!a_set){
        std::cout << "a not set!" << std::endl;
        exit(EXIT_FAILURE);
    }


}
bool Parameters::MC_cycles_set=false;
bool Parameters::alpha_max_set=false;
bool Parameters::alpha_min_set=false;
bool Parameters::alpha_num_set=false;
bool Parameters::beta_set=false;
bool Parameters::omega_set=false;
bool Parameters::omega_z_set=false;
bool Parameters::a_set=false;

int Parameters::text = 0;
int Parameters::MC_cycles = 0;
double Parameters::alpha_min = 0;
double Parameters::alpha_max = 0;
int Parameters::alpha_num = 0;
double Parameters::beta = 0;
double Parameters::omega = 0;
double Parameters::omega_z = 0;
double Parameters::a = 0;




