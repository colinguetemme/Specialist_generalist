
#include "environment.h"

environment::environment(){
}
environment::~environment(){
}

void environment::initialise(parameters para){	// Need to add a parameter

	std::random_device rd; //!//
	std::mt19937 gen(rd());

	int tmax;
	if (para.env.stochastic == 0.0){
		// when d == 0, the environment is cyclic with period L 
		// so we just need the values in one cycle
		tmax = para.env.period;
		tmax = 10000; //!// TEST
	} else {
		// when d > 0, the simulation stop at 50000 generations Burger(2002)
		tmax = 10000;
	}

	vector<double> opti_G(tmax, 0);

	
	double mean_env = (para.env.min + para.env.max) / 2;
	double range_env = para.env.max - para.env.min;

	// a rng following a normal distribution
	std::normal_distribution<> norm(0, para.env.var);

	for (int t = 0; t < tmax; t++)
	{
		opti_G[t] = mean_env + para.env.amplitude * range_env * sin(2 * M_PI * t / para.env.period);
	}
    optimum = opti_G;
}

