#pragma once

// The subclass containing the environment parameters
class para_env{
    public:

    para_env();
    ~para_env();
    double amplitude; // The amplitude of the variation
    int period; // The period of one cycle in number of generations
    double stochastic; // The stochasticity factor w/ 0 a non stochastic environment
    double var = stochastic * amplitude; // Variance of the normal distribution in wich the optimum is sampled
};

// The classe gathering all the parameters for the model
class parameters{
    public:

    // Constructor / Destructor
    parameters(double, double, int);
    ~parameters();

    ///////////////////
	// THE VARIABLES //
    ///////////////////

    para_env env; // The subclass containing the environment parameters
};

