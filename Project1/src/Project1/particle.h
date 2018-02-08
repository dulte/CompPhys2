#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vec3/vec3.h"
#include <vector>
#include <math.h>


class Particle
{
public:
    Particle(std::vector<double> m_r, double m_alpha, double m_beta);

    void accept_step();
    void update_alpha(double);
    double get_length();

    //vec3 r;
    std::vector<double> r;

    double r_norm;
    double r_squared;

    //vec3 next_r;
    std::vector<double> next_r;

    double mass = 1;

    double alpha;
    double beta;
};

#endif // PARTICLE_H
