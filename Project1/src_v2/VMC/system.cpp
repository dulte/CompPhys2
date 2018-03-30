#include "system.h"


System::System()
{

    r = Eigen::MatrixXd(Parameters::dimension,Parameters::N);
    next_r = Eigen::MatrixXd(Parameters::dimension,Parameters::N);

    distance.resize(Parameters::N,Parameters::N);
    next_distance.resize(Parameters::N,Parameters::N);

    quantum_force_vector.resize(Parameters::dimension);
    quantum_force_vector_new.resize(Parameters::dimension);

    std::random_device rd;
    gen = std::mt19937_64(rd());
    distribution = std::normal_distribution<double>(0.0,1.0);

    h=1e-5;

    if(Parameters::a !=0){
           //System::compute_energy_numeric = &System::calculate_energy_interacting;
           System::wavefunction_function_pointer=&System::update_wavefunction_interacting;
           System::compute_local_energy=&System::get_local_energy_interacting;

       }
     else{
           //System::compute_energy_numeric = &System::calculate_energy_interacting;//FIX
           System::wavefunction_function_pointer=&System::update_wavefunction_noninteracting;
           System::compute_local_energy=&System::get_local_energy_noninteracting;

    }

    //Sets the seed of rand() to the current time
    srand(time(NULL));

}

void System::make_grid(double m_alpha){
    alpha = m_alpha;
    //Sets all positions to a random position [-1,1]
    if(a!=0){
        distribute_particles_interacting();
    }
    else{
        distribute_particles_noninteracting();
    }

    next_r = r;
    update();
    wavefunction_value=get_wavefunction();
}

void System::distribute_particles_interacting(){
    bool safe_distance = false;
    for(int i = 0;i<N;i++){
        safe_distance = false;
        while(!safe_distance){
            for(int j = 0;j<dimension;j++){
                if(D!=0){
                    r(j,i) = distribution(gen)*sqrt(dx);
                }else{
                    r(j,i) = 2*((double)rand()/RAND_MAX - 0.5);
                }
            }
            safe_distance = true;
            for(int n = 0;n<i;n++){
                if((r.col(n) - r.col(i)).norm() <= a){
                    safe_distance = false;
                    break;
                }
            }
        }
        std::cout << "Particle " << i << " placed" << std::endl;
    }
}
void System::distribute_particles_noninteracting(){
    for(int i = 0;i<N;i++){
        for(int j = 0;j<dimension;j++){
            if(D!=0){
                r(j,i) = distribution(gen)*sqrt(dx);
            }else{
                r(j,i) = 2*((double)rand()/RAND_MAX - 0.5);
            }
        }
    }
}


//Updates the distances between the particles
void System::update(){
    double temp_value = 0;
    for(int i = 0; i<N;i++){
        for(int j = 0;j<i;j++){
            temp_value = (r.col(i)- r.col(j)).norm();
            distance(i,j) = temp_value;
            distance(j,i) = temp_value;
        }
    }

    next_distance = distance;
    if(D!=0){
        quantum_force(0);
    }
    update_expectation();
}

void System::make_move_and_update(const int move){
    //Makes a random move
    if(D!=0){
        quantum_force(move);
    }
    double random_nr = 0;
    for(int i = 0; i<dimension; i++){
        if(D!=0){
            random_nr = sqrt(dx)*distribution(gen) + quantum_force_vector(i)*dx*D;
        }else{
            random_nr = dx*2*((double)rand()/RAND_MAX - 0.5);
        }

        next_r(i,move) = r(i,move) +  random_nr;
    }
    update_next_distance(move);



}

void System::update_next_distance(int move){
    double dist = 0;
    for(int i = 0;i<N;i++){
        dist = (next_r.col(i)- next_r.col(move)).norm();

        next_distance(i,move) = dist;
        next_distance(move,i) = dist;
    }

}

void System::update_expectation(){
    expectation_derivative=0;
    expectation_derivative_energy=0;
    expectation_local_energy=0;
    expectation_local_energy_squared=0;
    wavefunction_probability=0;
    wavefunction_value=get_wavefunction();
}

double System::check_acceptance_and_return_energy(int move){
    //Random value [0,1]
    double temp_value = (double)rand()/RAND_MAX;

    //If r is less than the acceptance prob, r is updated to the new r
    if(temp_value <= get_probability_ratio(move)){
        update_wavefunction(move);
        r.col(move) = next_r.col(move);

        distance.col(move) = next_distance.col(move);
        distance.row(move) = next_distance.row(move);

    }
    else{
        next_r.col(move) = r.col(move);
        next_distance.col(move) = distance.col(move);
        next_distance.row(move) = distance.row(move);

    }
    //return get_local_energy_noninteracting();
    //return calculate_energy_interacting();
    //return get_local_energy();
    return 0;
}


double System::phi_exponant(const Eigen::VectorXd &r){
    temp_value = 0;

    for(int i = 0;i<dimension;i++){
        if(i == 2){
            //Multiplices beta to the z-componant
            temp_value += beta*r(i)*r(i);
        }
        else{
            temp_value += r(i)*r(i);
        }
    }
    return -alpha*temp_value;
}

double System::f(double dist){
    double function;
    if(dist <= a){
        function = 0;
    }
    else{
        function = 1 - a/dist;
    }

    return function;
}

double System::get_probability_ratio(int move){
    double temp_value = phi_exponant(r.col(move)); //Stores the probability before move
    double temp_value2 = phi_exponant(next_r.col(move)); //Stores the probability of move
    double f_part = 1;
    double green_part = 1;
    if(a!=0){
        f_part *= update_wavefunction_interacting_f(move);
    }
    if(D!=0){
        green_part = greens_function_ratio(move);
    }


    return exp(2*(temp_value2-temp_value))*f_part*f_part*green_part;
}

double System::get_wavefunction(){
    double temp_value = 0; //Stores the exponants of phi
    double f_part = 1;
    Eigen::VectorXd temp_r;
    for(int i = 0;i<N;i++){
        temp_r = r.col(i);
        temp_value += phi_exponant(temp_r);
        if(a != 0){
            for(int other = i+1;other<N;other++){
                f_part*= f(distance(other,i));
            }
        }

    }
    return exp(temp_value)*f_part;
}

void System::update_wavefunction(const int move){
    return (this->*wavefunction_function_pointer)(move);
}

void System::update_wavefunction_noninteracting(const int move){
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)));

}

void System::update_wavefunction_interacting(const int move){
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)))*update_wavefunction_interacting_f(move);

}


double System::update_wavefunction_interacting_f(const int move){
    double second_factor_of_psi=1;
    for(int i = 0; i < N; i++){
            if(i != move){
                 second_factor_of_psi*= f(next_distance(i,move))/f(distance(i,move));
        }
    }
    return second_factor_of_psi;
}

double System::get_probability(){
    double temp_value = get_wavefunction();
    return temp_value*temp_value;
}


double System::get_local_energy(){
    return (this->*compute_local_energy)();
}


double System::get_local_energy_noninteracting(){
        double total_energy = 0;
        double temp_value = 0;
        //temp_r = Eigen::VectorXd::Zero(dimension);
        //r_temp = Eigen::VectorXd::Zero(dimension);
        double factor1_noB = -2*(dimension)*alpha*N ;
        double factor1_B = -2*alpha*(dimension - 1)*N -  2*alpha*beta*N;
        double factor2 = 4*alpha*alpha;
        double pot_factor = 0.5*omega*omega;
        double r_i_annen = 0;
        double wavefunction_derivative_value=0;
        double omega_ratio = omega_z/omega;


        for(int k = 0;k<N;k++){
            for(int i = 0; i<dimension;i++){
                if(i==2){
                    temp_value += r(i,k)*r(i,k)*beta*beta;
                    wavefunction_derivative_value+=beta*r(i,k)*r(i,k);
                    r_i_annen += r(i,k)*r(i,k);//omega_ratio*r(i,k)*r(i,k);
                }
                else{
                    temp_value += r(i,k)*r(i,k);
                    wavefunction_derivative_value+=r(i,k)*r(i,k );
                    r_i_annen += r(i,k)*r(i,k);
                }


            }
        }

        if(dimension >= 3){
            total_energy = factor1_B + factor2*temp_value;
        }
        else{
            total_energy = factor1_noB + factor2*temp_value;
        }

        wavefunction_derivative_value*=-1;
        //temp_value = 0;
        temp_value=-0.5*total_energy+ pot_factor*r_i_annen;
        expectation_local_energy+=temp_value;
        expectation_local_energy_squared+=temp_value*temp_value;
        expectation_derivative+=wavefunction_derivative_value;
        expectation_derivative_energy+=(wavefunction_derivative_value)*temp_value;
        //std::cout<<temp_value<<std::endl;




    return temp_value;
}

double System::udiv(int idx1,int idx2){
    double du_dphi;
    if(distance(idx1,idx2) <= a){
        du_dphi = 1e15;
        std::cout << "Help it happend!!!" << std::endl;
        exit(1);
    }
    if(distance(idx1,idx2) > a){
        du_dphi = a/(distance(idx1,idx2)*distance(idx1,idx2) - a*distance(idx1,idx2));

    }
    return du_dphi;
    }

double System::udivdiv(int idx1,int idx2){
    double du2_dphi2;
    double temp_r = distance(idx1,idx2);
    if(distance(idx1,idx2) <= a){
        du2_dphi2 = 1e15;

    }
    if(distance(idx1,idx2) > a){
        du2_dphi2 = (a*a - 2*a*distance(idx1,idx2))/(
                    (distance(idx1,idx2)*distance(idx1,idx2)-
                     a*distance(idx1,idx2))*(distance(idx1,idx2)*distance(idx1,idx2)-a*distance(idx1,idx2)));
        //du2_dphi2 = -(a/((temp_r-a)*(temp_r-a))*a/(temp_r*temp_r) + 2*a/(temp_r*temp_r*temp_r - a*temp_r*temp_r));
    }
    return du2_dphi2;
}


double System::get_local_energy_interacting(){
    double sec_fac = 0;
    double trd_fac = 0;
    double frt_fac = 0;

    double total_energy = 0;
    double temp_value = 0;
    temp_r = Eigen::VectorXd::Zero(dimension);
    Eigen::VectorXd grad_r = Eigen::VectorXd::Zero(dimension);
    Eigen::VectorXd temp_difference_vector = Eigen::VectorXd::Zero(dimension);
    double factor1_noB = -2*(dimension)*alpha*N ;
    double factor1_B = -2*alpha*(dimension -1)*N -2*alpha*beta*N;
    double factor2 = 4*alpha*alpha;
    double pot_factor = 0.5*omega*omega;
    double r_i_annen = 0;
    double wavefunction_derivative_value=0;
    double omega_ratio = omega_z/omega;



    for(int k = 0;k<N;k++){
        for(int i = 0; i<dimension;i++){
            if(i==2){
                temp_value += r(i,k)*r(i,k)*beta*beta;
                r_i_annen += omega_ratio*r(i,k)*r(i,k);
                wavefunction_derivative_value+=beta*r(i,k)*r(i,k);
            }
            else{
                temp_value += r(i,k)*r(i,k);
                wavefunction_derivative_value+=r(i,k)*r(i,k);
                r_i_annen += r(i,k)*r(i,k);
            }
        }
    }


    for(int idx1 = 0; idx1 < N; idx1++ ){
        for(int idx2 = 0; idx2 < N; idx2++){
            for(int dim = 0; dim < dimension; dim++){
                if(dim == 2){

                    if(idx1 != idx2){
                        sec_fac += 2*udiv(idx1,idx2)*beta*beta*(((r(dim,idx1)*r(dim,idx1) -r(dim,idx1)*r(dim,idx2))
                                                     /distance(idx1,idx2)));
                        frt_fac += udivdiv(idx1,idx2) + 2*udiv(idx1,idx2)/distance(idx1,idx2);

                    }
                }
                else{
                    if(idx1 != idx2){

                        sec_fac += 2*udiv(idx1,idx2)*(((r(dim,idx1)*r(dim,idx1) -r(dim,idx1)*r(dim,idx2))
                                                     /distance(idx1,idx2)));

                        frt_fac += udivdiv(idx1,idx2) + 2*udiv(idx1,idx2)/distance(idx1,idx2);

                    }

                }
            }
        }
    }


    for(int idx11 = 0; idx11 < N; idx11++){
        for(int idx22 = 0; idx22 < N; idx22++){
            for(int idx33 = 0; idx33 < N; idx33++){
                for(int d = 0; d < dimension; d++){
                    if( idx11 != idx22 && idx11 !=idx33){
                        trd_fac += udiv(idx11,idx33)*udiv(idx11,idx22)*(r(d,idx11)*r(d,idx11)
                                                                        - r(d,idx11)*r(d,idx22)
                                                                        -r(d,idx11)*r(d,idx33)
                                                                        +r(d,idx22)*r(d,idx33))/(distance(idx11,idx22)*distance(idx11,idx33));
                    }
                }

            }
        }
    }


    /*

    for(int k = 0;k<N;k++){
        temp_r = r.col(k);
        for(int dim = 0;dim<dimension;dim++){
            if(dim == 2){
                grad_r(dim) = temp_r(dim)*beta;
            }
            else{
                grad_r(dim) = temp_r(dim);
            }
        }
        grad_r *= -4*alpha;
        for(int j = 0;j<N;j++){
            if(j!=k){
                temp_difference_vector += (r.col(k)-r.col(j))/(distance(k,j))*udiv(k,j);
            }
        }

        sec_fac += grad_r.dot(temp_difference_vector);
    }

    for(int k = 0; k<N;k++){
        for(int j = 0; j<N;j++){
            for(int i = 0; i<N;i++){
                if(k!=i && k!=j){
                    trd_fac += (r.col(k) - r.col(i)).dot(r.col(k) - r.col(j))/(distance(k,j)*distance(k,i))*udiv(k,i)*udiv(k,j);
                }
           }
        }
    }
    for(int k = 0;k<N;k++){
        for(int j = 0;j<N;j++){
            if(j!=k){
                frt_fac += udivdiv(k,j) + 2.0/(distance(k,j))*udiv(k,j);
            }
        }
    }*/


    sec_fac = -2*alpha*sec_fac;
    if(dimension >= 3){
        total_energy = factor1_B + factor2*temp_value + sec_fac + trd_fac + frt_fac;
    }
    else{
        total_energy = factor1_noB + factor2*temp_value + sec_fac + trd_fac + frt_fac;
    }

    temp_value= -0.5*total_energy+ pot_factor*r_i_annen;
    expectation_local_energy+=temp_value;
    expectation_local_energy_squared+=temp_value*temp_value;
    expectation_derivative+=wavefunction_derivative_value;
    expectation_derivative_energy+=(wavefunction_derivative_value)*temp_value;

    return temp_value;
}

void System::quantum_force(int move){
    double grad_value = 0;
    double grad_value_new = 0;




    for(int j=0; j<dimension;j++){

        if(a!=0){
            for(int k=0; k<N;k++){
                if(k !=move){
                    if(distance(move,k)>a){
                        grad_value+=2*a*(r(j,move)-(r(j,k)))/distance(move,k)*udiv(move,k);
                    }
                    if(next_distance(move,k)>a){
                        grad_value_new+=2*a*(next_r(j,move)-next_r(j,k))/next_distance(move,k)*udiv(move,k);
                    }
                 }
            }
        }
        if(j==2){
           quantum_force_vector(j)=-4.0*alpha*beta*r(j,move)+grad_value;
           quantum_force_vector_new(j)=-4.0*alpha*beta*next_r(j,move)+grad_value;

        }
        else{
            quantum_force_vector(j)=-4.0*alpha*r(j,move)+grad_value;
            quantum_force_vector_new(j)=-4.0*alpha*next_r(j,move)+grad_value;
        }
        grad_value=0;
        grad_value_new=0;
    }


}


double System::greens_function_ratio(int move)
{
    double value=0;
    double value_new=0;
    double exponent = 0;
    double exponent_new = 0;
    double exponent_factor=D*dx;

    //temp_r = Eigen::VectorXd::Zero(dimension);

    quantum_force(move);
    /*
    for(int i=0; i<N; i++){
        if(i==move){
            temp_r = next_r.col(move);
         }
         else{
             temp_r = r.col(move);
         }
         for(int k=0;k<dimension;k++){
                exponent += ((r(k,i)-temp_r(k))-exponent_factor*quantum_force_matrix(i,k))*((r(k,i)-temp_r(k))-exponent_factor*quantum_force_matrix(i,k));
                exponent_new += ((temp_r(k)-r(k,i))-exponent_factor*quantum_force_matrix_new(i,k))*((temp_r(k)-r(k,i))-exponent_factor*quantum_force_matrix_new(i,k));

         }

        exponent=exponent/(4.0*exponent_factor);
        exponent_new=exponent_new/(4.0*exponent_factor);
        value+=exp(-exponent);
        value_new+=exp(-exponent_new);

    }*/
/*
    for(int k=0;k<dimension;k++){
        value += 0.5*(quantum_force_vector(k) - quantum_force_vector_new(k))*(D*dx*0.5*(quantum_force_vector(k)-quantum_force_vector_new(k)) - next_r(k,move) + r(k,move));

    }*/

    value = -(r.col(move) - next_r.col(move) -D*dx*quantum_force_vector_new).squaredNorm()+(next_r.col(move)-r.col(move) - D*dx*quantum_force_vector).squaredNorm();
    value /= 4*D*dx;

    return exp(value)+N-1;
    //return value_new/value;

}


double System::calculate_energy_interacting(){
    double wavefunction_value_plus = 0;
    double wavefunction_value_minus = 0;

    double kinetic_energy = 0;
    double potential_energy = 0;
    wavefunction_value=get_wavefunction();
    for(int i = 0; i<N;i++){
        potential_energy+=omega*omega*r.col(i).squaredNorm();
        for(int j = 0; j<dimension;j++){
            //std::cout<<r->coeffRef(j,i)<<std::endl;
            r(j,i)+=h;
            wavefunction_value_plus=get_wavefunction();
            r(j,i)-=2*h;
            wavefunction_value_minus=get_wavefunction();
            r(j,i)+=h;
            kinetic_energy -= (wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value)/h/h/wavefunction_value;
            //std::cout<<r->coeffRef(j,i)<<std::endl;
        }
        //kinetic_energy += 2*dimension*alpha-4*alpha*alpha*r->col(i).squaredNorm();

    }
    //std::cout<<wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value<<std::endl;

    return (0.5*potential_energy+0.5*(kinetic_energy));
}





