#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace Parameters
{

    void read_parameters(std::string location);
    extern int text;
    extern int MC_cycles;


    extern double alpha_max;
    extern double alpha_min;
    extern int alpha_num;

    extern double beta;

    extern double omega;
}

#endif // PARAMETERS_H
