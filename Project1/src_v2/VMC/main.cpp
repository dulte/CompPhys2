#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "mpi.h"
#include "Parameters/parameters.h"
#include "simulation.h"
#include "system.h"
#include <string>

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
    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    Eigen::initParallel();
    //Initiates system and simulation
    System * system = new System();
    Simulation * simulation = new Simulation(system);

    std::string filename = "..//output//data";
    filename.append(std::to_string(my_rank));
    filename.append(".txt");

    simulation->initiate();
    //simulation->run();
    //std::cout << "Correct a: " << simulation->conjugate_gradient(0.5, 0.01)<<std::endl;
    simulation->run(filename);

    delete simulation;
    delete system;

    double EndTime = MPI_Wtime();
    double TotalTime = EndTime-StartTime;

    if ( my_rank == 0 )  cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;

    MPI_Finalize ();
    return 0;
}
