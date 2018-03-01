#include "system.h"


System::System()
{

    r = Eigen::MatrixXd(Parameters::dimension,Parameters::N);
    std::cout << "Ehi" << std::endl;
    next_r = Eigen::MatrixXd(Parameters::dimension,Parameters::N);
    distance.resize(Parameters::N,Parameters::N);
    next_distance.resize(Parameters::N,Parameters::N);
    quantum_force_matrix.resize(Parameters::N,Parameters::N);
    kinetic_energy=0;
    potential_energy=0;
    wavefunction_probability=0;
    wavefunction_value=get_wavefunction();
    h=1e-5;


    if(Parameters::a !=0){
        System::make_move = &System::make_move_and_update_interacting;
        System::compute_energy = &System::calculate_energy_interacting; //FIX
    }
    else{
        System::make_move=&System::make_move_and_update_non_interacting;
        System::compute_energy = &System::calculate_energy_interacting;//FIX
    }

    //Temp variables
    temp_r.resize(Parameters::dimension);
    temp_r2.resize(Parameters::dimension);
    temp_N_vec.resize(N);
    //temp_value = 0;
    //temp_value2 = 0;

    //Sets the seed of rand() to the current time
    srand(time(NULL));
}


void System::make_grid(double m_alpha){
    alpha = m_alpha;
    //Sets all positions to a random position [-1,1]
    //r = Eigen::MatrixXd::Random(Parameters::dimension,Parameters::N);


    r = Eigen::MatrixXd::Constant(Parameters::dimension,Parameters::N,0.0);
    next_r = r;
    update();
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
}


void System::make_move_and_update(const int move){
    (this->*make_move)(move);
}

double System::calculate_energy(){
    return (this->*compute_energy)();

}


void System::make_move_and_update_non_interacting(const int move){
    //Makes a random move
    for(int i = 0; i<dimension; i++){
        next_r(i,move) += dx*((double)rand()/RAND_MAX - 0.5);
    }
    std::cout << next_r << std::endl;
    std::cout << r << std::endl;
    //double temp_value = 0;
    //Updates the distance matrix after move
    for(int i = 0;i<N;i++){
        if(i != move){
            temp_value = (r->col(move)- r->col(i)).norm();
            next_distance(i,move) = temp_value;
            next_distance(move,i) = temp_value;
        }
    }
}

double System::check_acceptance_and_return_energy(int move){
    //Random value [0,1]
    double temp_value = (double)rand()/RAND_MAX;

    //If r is less than the acceptance prob, r is updated to the new r
    if(temp_value <= get_probability_ratio(move)){
        update_wavefunction(move);
        r = next_r;
        acceptance++;
        //update();
        //distance = next_distance;
    }
    else{
        next_r=r;
    }
    return get_local_energy();
}


double System::phi_exponant(const Eigen::VectorXd &r){
    double temp_value = 0;

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

    return exp(2*(temp_value2-temp_value));
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

double System::get_probability(){
    double temp_value = get_wavefunction();
    return temp_value*temp_value;
}

double System::calculate_energy_noninteracting(){
    kinetic_energy=0;
    potential_energy=0;
    for(int i=0; i<N;i++){
        temp_r=r->col(i);
        temp_r2=r->col(i);
        potential_energy += 0.5*omega*temp_r.squaredNorm();
        for(int j=0; j<dimension;j++){
            //temp_r2(j)+=h;
            //kinetic_energy+=exp(phi_exponant(temp_r2-temp_r));
            //temp_r2(j)-=2*h;
            //kinetic_energy+=exp(phi_exponant(temp_r2-temp_r));
            kinetic_energy+=exp(-2*alpha*h*temp_r(j)-alpha*h*h)+exp(2*alpha*temp_r(j)-alpha*h*h)-2;
        }
    }
    kinetic_energy *= -1/(2.0*h*h);
    std::cout<<kinetic_energy<<std::endl;
    return kinetic_energy+potential_energy;

}

double System::calculate_energy_interacting(){
    wavefunction_value_plus = 0;
    wavefunction_value_minus = 0;

    kinetic_energy = 0;
    potential_energy = 0;
    wavefunction_value=get_wavefunction();
    for(int i = 0; i<N;i++){
        potential_energy+=omega*omega*r->col(i).squaredNorm();
        for(int j = 0; j<dimension;j++){
            //std::cout<<r->coeffRef(j,i)<<std::endl;
            r->coeffRef(j,i)+=h;
            wavefunction_value_plus=get_wavefunction();
            r->coeffRef(j,i)-=2*h;
            wavefunction_value_minus=get_wavefunction();
            r->coeffRef(j,i)+=h;
            kinetic_energy -= (wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value)/h/h/wavefunction_value;
            //std::cout<<r->coeffRef(j,i)<<std::endl;
        }
        //kinetic_energy += 2*dimension*alpha-4*alpha*alpha*r->col(i).squaredNorm();

    }
    //sstd::cout<<wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value<<std::endl;

    return (0.5*potential_energy+0.5*(kinetic_energy));
}


void System::update_wavefunction(const int move){
    wavefunction_value*=exp(phi_exponant(next_r->col(move))-phi_exponant(r->col(move)));

}

double System::get_local_energy(){
    double total_energy = 0;
    double temp_value = 0;


    for(int k = 0;k<N;k++){
        for(int i = 0; i<dimension;i++){
            if(i==2){
                temp_value += r(i,k)*r(i,k)*beta*beta;
            }
            else{
                temp_value += r(i,k)*r(i,k);
            }
        }

        temp_value *= 4*alpha*alpha;
        temp_value -= 2*dimension*alpha;
        if(dimension>=3){
            temp_value += 2*alpha - 2*alpha*beta;
        }

        total_energy += temp_value;
        temp_value = 0;
    }

    return -0.5*total_energy;


    /*
    temp_value = new Eigen::VectorXd(dimension);

    for(int k = 0; k<N;k++){

        for(int i = 0;i<dimension;i++){
            if(i==2){
                temp_value += r->coeffRef(i,k)*r->coeffRef(i,k)*beta;
            }
            else{
                temp_value += r->coeffRef(i,k)*r->coeffRef(i,k);
            }
            temp_value *= -4*alpha*temp_value;
            temp_value -= 2*alpha;
            total_energy += temp_value;
            temp_value = 0;
        }

        for(int j = 0; j<N;j++){
            if(j!=k){
                temp_r += (r->col(k)-r->col(j))/distance(k,j);
            }

        }


    }
    */
    return 0;
}


double System::quantum_force(){
    double grad_value = 0;
    double temp_value = 0;
    for(int i=0; i<N;i++){
        for(int j=0; j<dimension;j++){
            for(int k=0; k<N;k++){
                if(k !=i){
                    if(distance(i,k)>a){
                        grad_value+=a*(r(j,i)-(r(j,k)))/(distance(i,k)*distance(i,k)*distance(i,k));
                    }
                 }
            }
            if(j==2){
               quantum_force_matrix(i,j)=-4.0*alpha*beta*r(j,i)+grad_value;

            }
            else{
                quantum_force_matrix(i,j)=-4.0*alpha*r(j,i)+grad_value;
            }
            grad_value=0;
        }
    }

}





