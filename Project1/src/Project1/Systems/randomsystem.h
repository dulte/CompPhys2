#ifndef RANDOMSYSTEM_H
#define RANDOMSYSTEM_H

#include "system.h"
#include <memory>
#include "Parameters/parameters.h"

class RandomSystem : public System
{
public:
    RandomSystem(std::shared_ptr<TrialFunction> m_trial);

    void update_alpha(double) override;
    void grid_setup(int,double) override;
    void propose_step() override;
    double check_acceptance_and_return_energy() override;

private:
    double step_size = 0.1;
    double beta = 1;
};

#endif // RANDOMSYSTEM_H
