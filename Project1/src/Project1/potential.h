#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Vec3/vec3.h"
#include <vector>
#include <iostream>

class Potential
{
public:
    //Potential();
    virtual double get_external_potential(double){};
    virtual double get_external_potential(std::vector<double>){};
    virtual double get_inter_potential(std::vector<double>,std::vector<double>){};

};

#endif // POTENTIAL_H
