#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "Parameters/parameters.h"
#include "simulation.h"
#include "system.h"
#include "mpi.h"

using namespace std;

int main(int nargs, char *args[])
{
    int numprocs, my_rank;
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    double StartTime = MPI_Wtime();

    if(my_rank == 0){

        cout << "MC steps, particles: " << Parameters::MC_cycles << ", " << Parameters::N << endl;
    }

    Eigen::initParallel();
    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    cout << "MC steps, particles: " << Parameters::MC_cycles << ", " << Parameters::N << endl;

    //Initiates system and simulation
    System * system = new System();
    Simulation * simulation = new Simulation(system);

    simulation->initiate();
    //double optimal_alpha = simulation->conjugate_gradient(0.3, 1.);
    //std::cout << "Correct a: " << optimal_alpha <<std::endl;
    simulation->run(my_rank);
    //simulation->data_for_derivated();
    //simulation->oneBodyDensity(optimal_alpha,0.03,0.,6.);


    double EndTime = MPI_Wtime();
    double TotalTime = EndTime-StartTime;

    if ( my_rank == 0 )  cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;

    MPI_Finalize ();
    return 0;
}
