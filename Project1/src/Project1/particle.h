#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vec3/vec3.h"


class Particle
{
public:
    Particle(vec3 m_r,double m_alpha,double m_beta);

    void accept_step();
    void update_alpha(double);

    vec3 r;
    double r_norm;
    double r_squared;

    vec3 next_r;

    double mass = 1;

    double alpha;
    double beta;
};

#endif // PARTICLE_H
