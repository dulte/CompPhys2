#include <iostream>
#include <memory>
#include "simulation.h"
#include "system.h"
#include "Systems/randomsystem.h"
#include "trialfunction.h"
#include "potential.h"
#include "Potentials/simpleharmonicoscillator.h"
#include "Parameters/parameters.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;

    //Read Parameters
    Parameters::read_parameters("../input/parameters.txt");
    cout << Parameters::text << endl;

    //Place all systems, potentials and trialfunctions here
    SimpleHarmonicOscillator potential = SimpleHarmonicOscillator(1);
    TrialFunction trial_function = TrialFunction(std::make_shared<SimpleHarmonicOscillator>(potential));
    RandomSystem system = RandomSystem(std::make_shared<TrialFunction>(trial_function));
    std::shared_ptr<System> pSystem = std::make_shared<RandomSystem>(system);


    Simulation simulation = Simulation(pSystem);
    simulation.initiate(1,0,1,1);

    simulation.run(1);

    cout << "Done!" << endl;

    return 0;
}
