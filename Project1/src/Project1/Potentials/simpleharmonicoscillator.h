#ifndef SIMPLEHARMONICOSCILLATOR_H
#define SIMPLEHARMONICOSCILLATOR_H

#include "../potential.h"
#include "../Parameters/parameters.h"

class SimpleHarmonicOscillator : public Potential
{
public:
    SimpleHarmonicOscillator(double);
    double get_external_potential(vec3);
    double get_inter_potential(vec3,vec3);

private:
    double omega = 1;
};

#endif // SIMPLEHARMONICOSCILLATOR_H
