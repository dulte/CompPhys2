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

    //Temp variables
    //temp_r.resize(Parameters::dimension);
    //temp_value = 0;
    //temp_value2 = 0;

    //Sets the seed of rand() to the current time
    srand(time(NULL));
}

void System::make_grid(double m_alpha){
    alpha = m_alpha;
        //Sets all positions to a random position [-1,1]
        //r = Eigen::MatrixXd::Random(Parameters::dimension,Parameters::N)*0.1;


        for(int i = 0;i<N;i++){
            for(int j = 0;j<dimension;j++){
                r(j,i) = distribution(gen)*sqrt(dx);
            }
        }

        //r = Eigen::MatrixXd::Constant(Parameters::dimension,Parameters::N,0.0);
    next_r = r;
    update();
    wavefunction_value=get_wavefunction();
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
    quantum_force(0);
}

void System::make_move_and_update(const int move){
    //Makes a random move
    double random_nr = 0;
    for(int i = 0; i<dimension; i++){
        random_nr = dx*distribution(gen);// + quantum_force_matrix(move,i)*dx*D;
        next_r(i,move) = r(i,move) +  random_nr;//((double)rand()/RAND_MAX - 0.5);
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

double System::check_acceptance_and_return_energy(int move){
    //Random value [0,1]
    double temp_value = (double)rand()/RAND_MAX;

    //If r is less than the acceptance prob, r is updated to the new r
    if(temp_value <= get_probability_ratio(move)){
        update_wavefunction(move);
        r = next_r;
        //update();
        distance = next_distance;
    }
    else{
        next_r = r;
    }
    return get_local_energy();
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

double System::get_probability_ratio(int move){
    double temp_value = phi_exponant(r.col(move)); //Stores the probability before move
    double temp_value2 = phi_exponant(next_r.col(move)); //Stores the probability of move

    return exp(2*(temp_value2-temp_value))*greens_function_ratio(move);
}

double System::get_wavefunction(){
    double temp_value = 0; //Stores the exponants of phi
    Eigen::VectorXd temp_r;
    for(int i = 0;i<N;i++){
        temp_r = r.col(i);
        temp_value += phi_exponant(temp_r);
    }
    return exp(temp_value);
}

void System::update_wavefunction(const int move){
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)));

}

double System::get_probability(){
    double temp_value = get_wavefunction();
    return temp_value*temp_value;
}


double System::get_local_energy(){
        double total_energy = 0;
        double temp_value = 0;
        //temp_r = Eigen::VectorXd::Zero(dimension);
        //r_temp = Eigen::VectorXd::Zero(dimension);
        double factor1_noB = -2*(dimension)*alpha*N ;
        double factor1_B = -2*alpha*(dimension - 1)*N -  2*alpha*beta*N;
        double factor2 = 4*alpha*alpha;
        double pot_factor = 0.5*omega*omega;
        double r_i_annen = 0;
        //double wavefunction_derivative_value=0;
        double omega_ratio = omega_z/omega;


        for(int k = 0;k<N;k++){
            for(int i = 0; i<dimension;i++){
                if(i==2){
                    temp_value += r(i,k)*r(i,k)*beta*beta;
                    //wavefunction_derivative_value+=beta*r(i,k)*r(i,k);
                    r_i_annen += r(i,k)*r(i,k);//omega_ratio*r(i,k)*r(i,k);
                }
                else{
                    temp_value += r(i,k)*r(i,k);
                    //wavefunction_derivative_value+=r(i,k)*r(i,k);
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

        //wavefunction_derivative_value*=-1;
        //temp_value = 0;
        temp_value=-0.5*total_energy+ pot_factor*r_i_annen;
        //expectation_local_energy+=temp_value;
        //expectation_local_energy_squared+=temp_value*temp_value;
        //expectation_derivative+=wavefunction_derivative_value;
        //expectation_derivative_energy+=(wavefunction_derivative_value)*temp_value;
        //std::cout<<temp_value<<std::endl;




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
                        grad_value+=2*a*(r(j,move)-(r(j,k)))/(distance(move,k)*distance(move,k)*distance(move,k));
                    }
                    if(next_distance(move,k)>a){
                        grad_value_new+=2*a*(next_r(j,move)-r(j,k))/(next_distance(move,k)*next_distance(move,k)*next_distance(move,k));
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

    temp_r = Eigen::VectorXd::Zero(dimension);

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

    for(int k=0;k<dimension;k++){
        value += 0.5*(quantum_force_vector(k) - quantum_force_vector_new(k))*(D*dx*0.5*(quantum_force_vector(k)-quantum_force_vector_new(k)) - next_r(k,move) + r(k,move));

    }
    return exp(value);
    //return value_new/value;

}





