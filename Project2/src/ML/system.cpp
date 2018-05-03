#include "system.h"

/*
This code runs the heavy-duty computational parts,
including updating the wavefunction and the energy.
*/
System::System()
{

    r = Eigen::MatrixXd(dimension,P);
    next_r = Eigen::MatrixXd(dimension,P);

    distance.resize(P,P);
    next_distance.resize(P,P);

    quantum_force_vector.resize(Parameters::dimension);
    quantum_force_vector_new.resize(Parameters::dimension);


    a_bias.resize(M);
    b_bias.resize(N);
    weights.resize(M,N);

    X.resize(M);
    X_next.resize(M);
    H.resize(N);

    std::random_device rd;
    gen = std::mt19937_64(rd());
    distribution = std::normal_distribution<double>(0.0,1.0);

    h=1e-6; //Step size for numerical derivative
    number_accept = 0;

    if(is_interacting){
	   //We use some function pointers to avoid some if statements in the functions
           System::wavefunction_function_pointer=&System::update_wavefunction_interacting;
           //System::compute_local_energy=&System::get_local_energy_interacting;

       }
     else{

           System::wavefunction_function_pointer=&System::update_wavefunction_noninteracting;
           System::compute_local_energy=&System::get_local_energy_noninteracting;

    }

    if(gibbs){
        System::gibbs_factor=0.5;
    }
    else{
        System::gibbs_factor=1.;
    }

    if(D!=0){
        System::greens_pointer=&System::greens_function_ratio;
    }
    else{
        System::greens_pointer=&System::greens_function_ratio_none;
    }

    //Sets the seed of rand() to the current time
    srand(time(NULL));
}

void System::make_grid(double m_alpha){
    /*
    This function initializes everything
    to prepare for a simulation with variational
    parameter m_alpha.
    */
    alpha = m_alpha;
    number_accept = 0;
    //Sets all positions to a random position [-1,1]

    distribute_particles_noninteracting();


    X_next = X;
    update();
    wavefunction_value=get_wavefunction();
}

void System::make_grid(Eigen::ArrayXd &parameters)
{

    Eigen::VectorXd par = (Eigen::VectorXd) parameters;
    Eigen::VectorXd w_flatten(M*N);
    for(int i = 0;i<par.size();i++){
        if(i<M){
            a_bias(i) = par(i);
        }else if(i>=M && i<(N+M)){
            b_bias(i-(M)) = par(i);
        }else{
            w_flatten(i-(M+N)) = par(i);
        }
    }

    Eigen::Map<Eigen::MatrixXd> M2(w_flatten.data(),M,N);

    weights = M2;



    distribute_particles_noninteracting();
    X_next = X;
    update();
    wavefunction_value=get_wavefunction();

}

void System::distribute_particles_interacting(){
    /*
    This function distributes the particles
    in the interacting case, ensuring that none
    of them start out closer than a from one another.
    */
    /*bool safe_distance = false;
    for(int i = 0;i<N;i++){
        safe_distance = false;
        while(!safe_distance){
            for(int j = 0;j<dimension;j++){
                if(D!=0){
                    r(j,i) = distribution(gen)*sqrt(dx);
                }else{
                    r(j,i) = 2*(static_cast<double>(rand())/RAND_MAX - 0.5);
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
    }*/
}
void System::distribute_particles_noninteracting(){
   /*
   This distribution is easier than the one above,
   as we do not need to take a minimum distance
   into account.
   */

    for(int i = 0;i<P*dimension;i++){
            if(D!=0){
                X(i) = distribution(gen)*sqrt(dx);
            }else{
                X(i) = 0.05*(static_cast<double>(rand())/RAND_MAX - 0.5);
        }
    }

}

void System::update_X_next(int move){
    double random_nr = 0;
    if(D!=0){
        quantum_force(move);
    }

    for(int i = 0; i<dimension; i++){
        if(D!=0){
            random_nr = sqrt(dx)*distribution(gen) + quantum_force_vector(i)*dx*D;
        }else{
            random_nr = dx*2*(static_cast<double>(rand())/RAND_MAX - 0.5);
        }

        X_next(move*dimension+i) = X(move*dimension+i) +  random_nr;
    }
 }



void System::update(){
    /*
    Updates the distance
    maxtrix after a new move
    */

    double dist = 0;
    for(int i = 0; i<P;i++){
        for(int j = 0;j<i;j++){
            for(int dim=0;dim<dimension;dim++){
                dist+=(X(i*dimension+dim)-X(j*dimension+dim))*(X(i*dimension+dim)-X(j*dimension+dim));
            }

            dist=sqrt(dist);
            distance(i,j) = dist;
            distance(j,i) = dist;
            dist = 0;
        }
        temp_value=0;
    }

    next_distance = distance;
    if(D!=0){
        quantum_force(0); //Re-initialize quantum force
    }
    update_expectation();
}

void System::make_move_and_update(const int move){
    /*
    Takes a random move as input and 
    proposes the new step.
    */

    update_X_next(move);
}

void System::update_next_distance(int move){
    /*
    Updates the next_distance matrix for
    to use in importance samplinf.
    */
    double dist = 0;
    for(int i = 0;i<P;i++){
        for(int dim=0;dim<dimension;dim++){
            dist+=(X_next(i*dimension+dim)-X_next(move*dimension+dim))*(X_next(i*dimension+dim)-X_next(move*dimension+dim));
        }
        dist=sqrt(dist);
        next_distance(i,move) = dist;
        next_distance(move,i) = dist;
        dist=0;
    }

}

void System::update_expectation(){
    /*
    Used for the gradient descent method
    */
    expectation_derivative=0;
    expectation_derivative_energy=0;
    expectation_local_energy=0;
    expectation_local_energy_squared=0;
    wavefunction_probability=0;
    wavefunction_value=get_wavefunction();
}

double System::check_acceptance_and_return_energy(int move){
    /*
    Performs the Metropolis test
    and updates the energy
    */


    //Random value [0,1]
    double temp_value = (double)rand()/RAND_MAX;

    //If r is less than the acceptance prob, r is updated to the new r

    if(temp_value <= get_probability_ratio(move)){

        //update_wavefunction(move);
        wavefunction_value = get_wavefunction();
        X = X_next;

        distance.col(move) = next_distance.col(move);
        distance.row(move) = next_distance.row(move);

        number_accept++;



    }
    else{

        X_next = X;
        next_distance.col(move) = distance.col(move);
        next_distance.row(move) = distance.row(move);


    }

    if(is_numerical){
        return calculate_energy_numerically();
    }
    else{
        //std::cout<<"Numerical energy: "<< calculate_energy_numerically()<<std::endl;
        return get_local_energy_noninteracting();
    }


}

double System::gibbs_sample_and_return_energy(){
    sample_x();
    sample_h();
    if(is_interacting){
        update();
    }
    return get_local_energy_noninteracting();
}


double System::phi_exponant(const Eigen::VectorXd &r){
    /*
    Smart way to update wavefunction
    */
    temp_value = 0;

    for(int i = 0;i<dimension;i++){

        temp_value += r(i)*r(i);

    }
    return -alpha*temp_value;
}



double System::f(double dist){
    /*
    Computes Jastrow factor
    */
    /*double function;
    if(dist <= a){
        function = 0;
    }
    else{
        function = 1 - a/dist;
    }*/

    return 1;//function;
}

double System::get_probability_ratio(int move){
    /*
    Smart way to compute probability ratio, 
    which accounts for interacting and importance sampling
    */

    double wavefunction_old=get_wavefunction();
    double wavefunction_new=get_wavefunction_next();

    return (wavefunction_new*wavefunction_new)/(wavefunction_old*wavefunction_old)*greens_function_ratio(move);

    /*
    double first_part_ratio = 0;
    double wave_function_second_part = 1;
    double wave_function_second_part_new = 1;
    double exp_factor = 0;
    double exp_factor_new = 0;

    first_part_ratio = -((X_next(move) - a_bias(move))*(X_next(move) - a_bias(move)) - (X(move) - a_bias(move))*(X(move) - a_bias(move)));
    first_part_ratio = exp(first_part_ratio/(2*sigma_squared));

    for(int j = 0;j<N;j++){
        exp_factor = 0;
        exp_factor_new = 0;
        for(int i=0;i<M;i++){
            exp_factor += X[i]*weights(i,j);
            exp_factor_new += X_next[i]*weights(i,j);

        }
        wave_function_second_part *= (1+exp(b_bias(j)+(1.0/sigma_squared)*exp_factor));
        wave_function_second_part_new *= (1+exp(b_bias(j)+(1.0/sigma_squared)*exp_factor_new));
    }

    return first_part_ratio*(wave_function_second_part_new/wave_function_second_part)*first_part_ratio*(wave_function_second_part_new/wave_function_second_part)*greens_factor(move);
    */

}





double System::get_wavefunction(){
    /*
    Computes the wavefunction (complete)
    */
    double wave_function_first_part = 0;
    double wave_function_second_part = 1;
    double exp_factor=0;

    for(int i=0;i<M;i++){
        wave_function_first_part+=-(X(i) - a_bias(i))*(X(i) - a_bias(i));
    }

    wave_function_first_part = exp(wave_function_first_part/(2*sigma_squared));

    for(int j = 0;j<N;j++){
        exp_factor=0;
        for(int i=0;i<M;i++){
            exp_factor+=X(i)*weights(i,j);
        }
        wave_function_second_part *= (1+exp(b_bias(j)+(1.0/sigma_squared)*exp_factor));
    }

    //std::cout << wave_function_first_part*wave_function_second_part << std::endl;

    if(gibbs){
        return sqrt(wave_function_first_part*wave_function_second_part);
    }
    else{
        return wave_function_first_part*wave_function_second_part;
    }
}


double System::get_wavefunction_next(){
    /*
    Computes the wavefunction (complete)
    */
    double wave_function_first_part = 0;
    double wave_function_second_part = 1;
    double exp_factor=0;

    for(int i=0;i<M;i++){
        wave_function_first_part+=-(X_next(i) - a_bias(i))*(X_next(i) - a_bias(i));
    }

    wave_function_first_part = exp(wave_function_first_part/(2*sigma_squared));

    for(int j = 0;j<N;j++){
        exp_factor=0;
        for(int i=0;i<M;i++){
            exp_factor+=X_next(i)*weights(i,j);
        }
        wave_function_second_part *= (1+exp(b_bias(j)+(1.0/sigma_squared)*exp_factor));
    }

    //std::cout << wave_function_first_part*wave_function_second_part << std::endl;

    if(gibbs){
        return sqrt(wave_function_first_part*wave_function_second_part);
    }
    else{
        return wave_function_first_part*wave_function_second_part;
    }
}




void System::update_wavefunction(const int move){
    return (this->*wavefunction_function_pointer)(move);
}

void System::update_wavefunction_noninteracting(const int move){
    /*
    Smart way to update wavefunction given move
    */
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)));

}

void System::update_wavefunction_interacting(const int move){
    /*
    Smart way to update wavefunction given move
    */
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)))*update_wavefunction_interacting_f(move);

}


double System::update_wavefunction_interacting_f(const int move){
    /*
    Returns the ratio of the jastrow Factor
    */
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

//Returns the local energy. Does also compute the derivative of E_L used in the gradien descent.
double System::get_local_energy_noninteracting(){
    double derivative_of_log_psi=0;
    double second_derivative_of_log_psi=0;
    double exp_factor=0;
    double denominator_factor=0;
    double potential_energy=0;
    double repulsive_interaction=0;
    double potential_factor=0.5*omega*omega;
    Eigen::VectorXd x_weight_product(N);



    for(int k=0;k<M;k++){
           derivative_of_log_psi+=(a_bias(k)-X(k))/(sigma_squared);
           second_derivative_of_log_psi+=-1.0/(sigma_squared);

           x_weight_product=(1.0/sigma_squared)*(weights.transpose()*X);

           for(int j=0;j<N;j++){
                exp_factor=exp(-b_bias(j)-x_weight_product(j));
                derivative_of_log_psi+=weights(k,j)/(sigma_squared*(1+exp_factor));
                denominator_factor = sigma_squared*sigma_squared*(1+exp_factor)*(1+exp_factor);
                second_derivative_of_log_psi+=(weights(k,j)*weights(k,j))*exp_factor/denominator_factor;
           }
        }

    if(is_interacting){
        for(int i=0;i<P;i++){
            for(int j=0;j<i;j++){
                repulsive_interaction+=1.0/distance(i,j);
            }
        }
    }

    for(int k=0;k<M;k+=dimension){
        for(int dim=0;dim<dimension;dim++)
            potential_energy+=X(k+dim)*X(k+dim);
    }

   return -0.5*(derivative_of_log_psi*derivative_of_log_psi*gibbs_factor*gibbs_factor+second_derivative_of_log_psi*gibbs_factor)+potential_factor*potential_energy+repulsive_interaction;
}

double System::d_psi_da(int k){
    return (X[k]-a_bias[k])/sigma_squared;
}

double System::d_psi_db(int k){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,k);
    }
    return gibbs_factor*1.0/(1+exp(-b_bias(k)-(1.0/sigma_squared)*exp_factor));
}

double System::d_psi_dw(int k, int l){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,l);
    }
    return gibbs_factor*X(k)/(sigma_squared*(1+exp(-b_bias(l)-(1.0/sigma_squared)*exp_factor)));
}

double System::d_psi_da_log(int k){
    return gibbs_factor*(X(k)-a_bias(k))/(2.0*sigma_squared);
}

double System::d_psi_db_log(int k){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,k);
    }
    return 1.0/(2.0*(1+exp(-b_bias(k)-(1.0/sigma_squared)*exp_factor)));
}

double System::d_psi_dw_log(int k, int l){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,l);
    }
    return gibbs_factor*X(k)/(2.0*sigma_squared*(1+exp(-b_bias(l)-(1.0/sigma_squared)*exp_factor)));
}


//Computes the quantum force
void System::quantum_force(int move){
    double grad_value = 0;
    double grad_value_next = 0;
    double exp_factor = 0;
    double exp_factor_next = 0;
    Eigen::VectorXd x_weight_product(N);
    Eigen::VectorXd x_weight_product_next(N);


    x_weight_product=(1.0/sigma_squared)*(weights.transpose()*X);
    x_weight_product_next=(1.0/sigma_squared)*(weights.transpose()*X_next);


    for(int dim=0;dim<dimension;dim++){
        grad_value+=(a_bias(move*dimension+dim)-X(move*dimension+dim));
        grad_value_next+=(a_bias(move*dimension+dim)-X_next(move*dimension+dim));
        grad_value/=sigma_squared;
        grad_value_next/=sigma_squared;
        for(int j=0;j<N;j++){
            exp_factor=exp(-b_bias(j)-x_weight_product(j));
            exp_factor_next=exp(-b_bias(j)-x_weight_product_next(j));
            grad_value+=weights(move,j)/(sigma_squared*(1+exp_factor));
            grad_value_next+=weights(move,j)/(sigma_squared*(1+exp_factor_next));
       }
       quantum_force_vector(dim)=gibbs_factor*grad_value;
       quantum_force_vector_new(dim)=gibbs_factor*grad_value_next;
       grad_value=0;
       grad_value_next=0;

    }
}

//Returns the ratio of the green functions
double System::greens_function_ratio(int move)
{
    double value=0;
    Eigen::VectorXd x_move(dimension);
    Eigen::VectorXd x_move_next(dimension);

    for(int dim=0;dim<dimension;dim++){
        x_move(dim)=X(dimension*move+dim);
        x_move_next(dim)=X_next(dimension*move+dim);

    }


    quantum_force(move);

    value = -(x_move - x_move_next -D*dx*quantum_force_vector_new).squaredNorm()+(x_move_next-x_move - D*dx*quantum_force_vector).squaredNorm();
    value /= 4*D*dx;

    return exp(value);
}

double System::greens_factor(const int move){
    return (this->*greens_pointer)(move);
}

double System::greens_function_ratio_none(int move){
    return 1.0;
}



//Calculates the energy numerically. Does NOT calculate the derivative of E_L
double System::calculate_energy_numerically(){
    double wavefunction_value_plus = 0;
    double wavefunction_value_minus = 0;

    double kinetic_energy = 0;
    double potential_energy = 0;
    double omega_ratio = 1;// omega_z/omega;
    wavefunction_value=get_wavefunction();
    for(int i = 0; i<M;i++){
            potential_energy+=omega*omega*X(i)*X(i);


            X(i)+=h;
            wavefunction_value_plus=get_wavefunction();
            X(i)-=2*h;
            wavefunction_value_minus=get_wavefunction();
            X(i)+=h;
            kinetic_energy -= (wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value)/h/h/wavefunction_value;

    }
    return (0.5*potential_energy+0.5*(kinetic_energy));
}



void System::sample_h(){
    std::uniform_real_distribution<double> sampling_distribution;

    Eigen::VectorXd backwards_H = -b_bias - (X.transpose()*weights).transpose()/sigma_squared;

    double coin_toss = 0;
    double sample_prob = 0;
    for(int i = 0;i<N;i++){
        sample_prob = 1./(1+exp(backwards_H(i)));
        coin_toss = 2*(static_cast<double>(rand())/RAND_MAX - 0.5);
        if(sample_prob > coin_toss){
            H(i) = 1.;
        }
        else{
            H(i) = 0.;
        }
    }
}

void System::sample_x(){
    std::normal_distribution<double> sampling_distribution;
    Eigen::VectorXd backwards_X = a_bias + weights*H;

    for(int i = 0;i<M;i++){
        sampling_distribution = std::normal_distribution<double>(backwards_X(i),sigma_squared);
        X(i) = sampling_distribution(gen);
    }
}




