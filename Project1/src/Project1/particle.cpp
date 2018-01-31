#include "particle.h"

Particle::Particle(std::vector<double>& m_r,double m_alpha,double m_beta)
{

    r = m_r;
    next_r = m_r;

    r_norm = get_length();
    r_squared = r_norm*r_norm;

    alpha = m_alpha;
    beta = m_beta;

}

void Particle::accept_step(){
    r = next_r;
    r_norm = get_length();
    r_squared = r_norm*r_norm;
}

void Particle::update_alpha(double m_alpha){
    alpha = m_alpha;
}

double Particle::get_length(){
    int size = r.size();
    double val = 0;

    for(int i = 0; i<size;i++){
        val = val + r[i]*r[i];
    }
    return sqrt(val);
}




