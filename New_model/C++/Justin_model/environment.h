#pragma once

# include <math.h>
# include <vector>
# include <random>

#include "parameters.h"

using namespace std;


// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());
std::uniform_real_distribution<> unif(0.0, 1.0);

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

    // Functions
    /**
     * @brief initialise the environment, maybe will be change to the constructor
     * @param param an object of type parameter as define in 'parameters.h'
     */
    void initialise(parameters);

    void output();
}chromo