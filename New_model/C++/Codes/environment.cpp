
#include "environment.h"

environment::environment(){
}
environment::~environment(){
}


/**
 * @brief Initialise the environment, 
 * 
 * @param para the 
 */
void environment::initialise(para_env para){	// Need to add a parameter

	std::random_device rd; //!//
	std::mt19937 gen(rd());

	int tmax;
	if (para.stochastic == 0.0){
		// when d == 0, the environment is cyclic with period L 
		// so we just need the values in one cycle and loop it
		tmax = 10000; //!// TEST
	} else {
		// when d > 0, the simulation stop at xx generations
		tmax = 10000;
	}

	vector<double> opti_G(tmax, 0);

	// For now since only using sin(), we ca describe 
	double mean_env = (para.min + para.max) / 2;
	double range_env = para.max - para.min;

	// a rng following a normal distribution
	std::normal_distribution<> norm(0, para.var);

	for (int t = 0; t < tmax; t++)
	{
		opti_G[t] = mean_env + para.amplitude * range_env * sin(2 * M_PI * t / para.period);
	}
    optimum = opti_G;
}

void environment::output(string out_file){
	std::ofstream outfile (out_file);
	for (int i = 0; i<optimum.size(); i++){
		outfile << i+1 << "," << optimum[i] << "\n";
	}
    outfile.close();
}

