#ifndef SIMPLEHARMONICOSCILLATOR_H
#define SIMPLEHARMONICOSCILLATOR_H

#include "../potential.h"
#include "../Parameters/parameters.h"
#include <vector>

class SimpleHarmonicOscillator : public Potential
{
public:
    SimpleHarmonicOscillator(double);
    double get_external_potential(double);
    double get_inter_potential(std::vector<double>,std::vector<double>);

private:
    double omega = 1;
};

#endif // SIMPLEHARMONICOSCILLATOR_H
