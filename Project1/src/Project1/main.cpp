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
    //Parameters params = Parameters::instance();
    Parameters::read_parameters("../input/parameters.txt");
    cout << Parameters::beta << endl;

    //Place all systems, potentials and trialfunctions here
    std::unique_ptr<HarmonicOscillator> potential = std::make_unique<HarmonicOscillator>();
    std::unique_ptr<TrialFunction> trial_function = std::make_unique<TrialFunction>(potential.get());
    std::unique_ptr<RandomSystem> system = std::make_unique<RandomSystem>(trial_function.get());
    //std::shared_ptr<System*> pSystem = std::make_shared<RandomSystem*>(system);


    std::shared_ptr<Simulation> simulation = std::make_shared<Simulation>(system.get());
    simulation->initiate();

    simulation->run(Parameters::MC_cycles);

    cout << "Done!" << endl;

    return 0;
}
