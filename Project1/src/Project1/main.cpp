#include <iostream>
#include <memory>
#include "simulation.h"
#include "system.h"
#include "Systems/randomsystem.h"
#include "trialfunction.h"
#include "potential.h"
#include "Potentials/harmonicoscillator.h"
#include "Parameters/parameters.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;

    //Read Parameters
    Parameters::read_parameters("../input/parameters.txt");

    //Place all systems, potentials and trialfunctions here
    HarmonicOscillator *potential = new HarmonicOscillator();
    TrialFunction * trial_function = new TrialFunction(potential);
    RandomSystem * system = new RandomSystem(trial_function);

    Simulation *simulation = new Simulation(system);
    simulation->initiate();

    simulation->run(Parameters::MC_cycles);

    cout << "Done!" << endl;

    delete simulation;
    delete system;
    delete trial_function;
    delete potential;

    return 0;
}
