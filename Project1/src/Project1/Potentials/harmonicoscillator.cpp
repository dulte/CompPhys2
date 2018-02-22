#include "harmonicoscillator.h"
#include "../potential.h"
#include "../Parameters/parameters.h"

HarmonicOscillator::HarmonicOscillator()
{
    omega = Parameters::omega;
    omega_z = Parameters::omega_z;
    dimension = Parameters::dimension;
    a = Parameters::a;
}

double HarmonicOscillator::get_external_potential(double m_r){
    return 0.5*omega*omega*m_r*m_r;
}

double HarmonicOscillator::get_external_potential(std::vector<double> r){
    double potential=0;

    for(int i=0;i<dimension;i++){

        if(i==2){
            potential+=omega_z*omega_z*(r[i]*r[i]);
        }
        else{
            potential+=omega*omega*(r[i]*r[i]);

        }
    }
    return 0.5*potential;

}

double HarmonicOscillator::get_inter_potential(std::vector<double> r_i, std::vector<double> r_j){
    double dist=0;
    for(int k=0;k<3;k++){
        dist+=(r_i[k]-r_j[k])*(r_i[k]-r_j[k]);
    }
    dist=sqrt(dist);
    if(dist<=a){
        return FLOAT_MAX;
    }
    else{
        return 0;
    }
}


