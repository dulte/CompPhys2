#include <iostream>
#include "simulation.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;

    Simulation simulation = Simulation();
    simulation.initiate(10,0,1,10);

    simulation.run(10);

    cout << "Done!" << endl;

    return 0;
}
