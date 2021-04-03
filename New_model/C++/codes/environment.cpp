
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
void environment::initialise_sinus(para_env para){	// Need to add a parameter

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
	double var = range_env*para.stochastic;
	std::normal_distribution<> norm(0, var);
	
	if (para.stochastic == 0){
		for (int t = 0; t < tmax; t++)
		{
			opti_G[t] = mean_env + para.amplitude * range_env * sin(2 * M_PI * t / para.period);
		}
	} else {
		for (int t = 0; t < tmax; t++)
		{
			opti_G[t] = mean_env + para.amplitude * range_env * sin(2 * M_PI * t / para.period) + norm(gen);
		}
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

/**
 * @brief More flexible environment, not perfect but gives the possibility to have non sinus function.
 * //!// The period should be changed to have a period rather than a time in each state.
 * For the stochasticity get 50000 steps and for the non-stochastic just generate one cycle (L*states.size()).
 * 
 * @param para, a para_env class containing the parameters of the environment
 * @param states vector containing the states of the changing environment
 *
 * @return The optimum genotype at each time step
 */
	void environment::initialise_custom(para_env para, vector<double> states){ 
	std::random_device rd; //!//
	std::mt19937 gen(rd());

	// The period should be a multiple of the number of steps
	if (para.period % states.size() != 0){
		para.period == states.size();
		cout << "The period is not a multiple of the number of states" << endl;
		cout << "The period has been change to " << states.size() << " (the number of states)" << endl;
		cout << " \n Continue by entering a number" << endl;
		int a = 0;
		cin >> a;
	}

	int L = para.period / states.size();
	vector<double> opti_G(para.tmax, 0);

	if (para.stochastic > 0){

		double mean_states = 0;
		double var_states = 0;
		for(int i = 0; i < states.size(); i++){
			mean_states += states[i];
		}
		mean_states /= states.size();
		for(int i = 0; i < states.size(); i++){
			var_states += pow((states[i]-mean_states),2);
		}
		double var_stoch = var_states * para.stochastic;
		std::normal_distribution<> norm(0, var_stoch);

		int index;
		for (int i = 0; i < para.tmax; i++){
			index = floor(i%(states.size()*L)/L); // The time step in a cycle
			opti_G[i] = states[index] + norm(gen);
		}

	} else {
		int index;
		for (int i = 0; i < para.tmax; i++){
			index = floor(i%(states.size()*L)/L); // The time step in a cycle
			opti_G[i] = states[index];
		}
	}
	optimum = opti_G;
}


