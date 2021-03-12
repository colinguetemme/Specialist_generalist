
#include "environment.h"
#include "parameters.h"

void environment::initialise(parameters param){	// Need to add a parameters
	int tmax;
	if (param.env.stochastic == 0.0){
		// when d == 0, the environment is cyclic with period L 
		// so we just need the values in one cycle
		tmax = param.env.period;
	} else {
		// when d > 0, the simulation stop at 50000 generations Burger(2002)
		tmax = 50000;
	}

	vector<double> opti_G(tmax, 0);

	// a rng following a normal distribution
	std::normal_distribution<> norm(0, param.env.var);

	for (int t = 0; t < tmax; t++)
	{
		opti_G[t] = norm(rdgen) + 0.5 + param.env.amplitude * sin(2 * M_PI * t / param.env.period);
	}
    optimum = opti_G;
}

