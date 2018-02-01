#include "harmonicoscillator.h"
#include "potential.h"
#include "../Parameters/parameters.h"

HarmonicOscillator::HarmonicOscillator()
{
    double omega = Parameters::omega;
    double omega_z = Parameters::omega_z;
    double a = Parameters::a;
}

double HarmonicOscillator::get_external_potential(double m_r){
    return 0.5*omega*omega*m_r*m_r;
}

double HarmonicOscillator::get_external_potential(std::vector<double> r){
    return 0.5*(omega*omega*(r.at(0)*r.at(0)+r.at(1)*r.at(0))+omega_z*omega_z*(r.at(2)+r.at(2)));

}

double HarmonicOscillator::get_inter_potential(std::vector<double> r_i, std::vector<double> r_j){
    double dist=0;
    for(int k=0;k<3;k++){
        dist+=(r_i[k]-r_j[k])*(r_i[k]-r_j[k]);
    }
    dist=sqrt(dist);
    if(dist<a){
        return FLOAT_MAX;
    }
    else{
        return 0;
    }
}


