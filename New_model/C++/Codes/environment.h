
#pragma once

# include <math.h>
# include <vector>
# include <random>

#include "parameters.h"
#include "model.h"

using namespace std;

// Random numbers generators
// seed random number generator

//std::uniform_real_distribution<> unif(0.0, 1.0);

/**
 * @brief The class that will store the environment
 * 
 */
class environment{
    public:

    // Constructor / Destructor
    environment();
    ~environment();

    // Variables
    vector<double> optimum; // Vector of optimum genotype fluctuating in time;

    // METHODS
    /**
     * @brief initialise the environment, maybe will be change to the constructor
     * 
     * Need to initialise the parameters before, parameter will include the amplitude, the period and
     * the stochasticity
     * 
     * @param param an object of type parameter as define in 'parameters.h'
     * 
     */
    void initialise(para_env);

    void output(string);
};