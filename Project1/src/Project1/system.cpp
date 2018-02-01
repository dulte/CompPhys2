#include "system.h"




System::System()
{
    size = Parameters::N;
    for(int i = 0;i<size;i++){
        rs.push_back(0);
    }
}



std::vector<double> System::get_postions(){
    for(int i = 0;i<size;i++){
        rs[i] = particles[i].r_norm;
    }

    return rs;
}
