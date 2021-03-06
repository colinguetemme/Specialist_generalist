#pragma once

# include <math.h>
# include <vector>
# include <random>

using namespace std;

// The structure containing the environment parameters
struct para_env{
    double min; // The minimum value of the environment
    double max; // The maximum value of the environment

    double amplitude; // The amplitude of the variation
    int period; // The period of one cycle in number of generations
    double stochastic; // The stochasticity factor w/ 0 a non stochastic environment
    double var = stochastic * amplitude; // Variance of the normal distribution in wich the optimum is sampled
};

// The structure containing the individual parameters
struct para_ind{

    double mut_rate; // The mutation rate for each individual
    // The pmf of the number of mutation (from mut rate) (poisson)
    std::poisson_distribution<> distr_mut_rat;

    double mean_val_mut; // The mean value of a mutation (should be negative)
    double var_val_mut; // The variance of the mutation value 

    // The pdf of the the fitness value of a mutation following a normal law (?)
    std::normal_distribution<> distr_mut_value; 
    
};

// The class gathering all the parameters for the model
class parameters{
    public:

    // Constructor / Destructor
    parameters(
        double, double, int, double, double, // environment
        double, double, double // individual
        );
    ~parameters();

    ///////////////////
	// THE VARIABLES //
    ///////////////////

    para_env env; // contains the parameters of the environment
    para_ind ind; // contains the parameters of the individuals (or pop)
};

