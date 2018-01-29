#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Parameters
{
public:
    static Parameters* m_instance;

    static void read_parameters(std::string location);
    static int text;
    static int MC_cycles;


    static double alpha_max;
    static double alpha_min;
    static int alpha_num;

    static double beta;

    static double omega;
};

#endif // PARAMETERS_H
