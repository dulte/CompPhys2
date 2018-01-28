#include "parameters.h"



int text;
int MC_cycles;
double alpha_max;
double alpha_min;
int alpha_num;

double beta;

double omega;

void Parameters::read_parameters(std::string location){
    std::ifstream infile;
    std::string line;

    infile.open(location);
    if(infile.fail()){
        std::cout << "Could not find parameter file at " << location << std::endl;
    }
    infile.clear();
    infile.seekg(0,std::ios::beg);
    for (std::string l; getline(infile,l);)
    {
        std::stringstream ss(l);

        std::string name;
        double variable;

        ss >> name >> variable;
        std::cout << name << std::endl;

        if (&name[0] == "#")
            continue;
        else if(name == "text"){
            text = variable;
        }
        else if(name == "MC_cycles"){
            MC_cycles = variable;
        }
        else if(name == "alpha_min"){
            alpha_min = variable;
        }
        else if(name == "alpha_max"){
            alpha_max = variable;
        }
        else if(name == "alpha_num"){
            alpha_num = variable;
        }

        else if(name == "beta"){
            beta = variable;
        }
        else if(name == "omega"){
            omega = variable;
        }
        else{
            std::cout << "Unknonw Variable found: " << name << std::endl;
        }

    }

}



