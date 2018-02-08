#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H

#include "../potential.h"
#include "../Parameters/parameters.h"
#include <vector>
#include <math.h>
#include <limits>


class HarmonicOscillator : public Potential
{
public:
    HarmonicOscillator();
    double get_external_potential(double);
    double get_external_potential(std::vector<double>);
    double get_inter_potential(std::vector<double>,std::vector<double>);
    double omega;
    double omega_z;
    double dimension;
    double a;

private:
    float FLOAT_MAX = std::numeric_limits<float>::max();

};

#endif // HARMONICOSCILLATOR_H
