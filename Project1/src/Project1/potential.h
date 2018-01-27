#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Vec3/vec3.h"

class Potential
{
public:
    //Potential();
    virtual double get_external_potential(vec3){};
    virtual double get_inter_potential(vec3,vec3){};

};

#endif // POTENTIAL_H
