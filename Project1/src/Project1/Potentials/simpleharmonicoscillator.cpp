#include "simpleharmonicoscillator.h"
#include "potential.h"


SimpleHarmonicOscillator::SimpleHarmonicOscillator(double m_omega)
{
    omega = m_omega;
}

double SimpleHarmonicOscillator::get_external_potential(vec3 m_r){
    return 0.5*omega*omega*m_r*m_r;
}

double SimpleHarmonicOscillator::get_inter_potential(vec3 r_i, vec3 r_j){
    return 0.0;
}


