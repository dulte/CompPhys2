#include "particle.h"

Particle::Particle(vec3 m_r,double m_alpha,double m_beta,double m_phi(vec3,double,double))
{

    r = m_r;
    next_r = m_r;

    phi = m_phi;

    alpha = m_alpha;
    beta = m_beta;
    phi_value = phi(r,alpha,beta);

}

void Particle::accept_step(){
    r = next_r;
    r_norm = r.length();
    r_squared = r_norm*r_norm;
    phi_value = phi(r,alpha,beta);
}

void Particle::update_alpha(double m_alpha){
    alpha = m_alpha;
}


