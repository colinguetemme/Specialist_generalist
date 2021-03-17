
# include "individual.h"

/**
 * @brief Calculate the fitness value of that individual from the optimum genotype of the environment
 * w
 * @param env_value The value of the tracked environmental parameter (in our experiment the temperature)
 */
void individual::fit_value(double env_value){
    fit_val = 1; // or init_fit_val if needed
    std::map<double, mutation>::iterator iter;
 
    for(iter = chr.mutations.begin(); iter != chr.mutations.end(); iter++){

        if(iter->second.min < env_value && iter->second.max > env_value){

           fit_val += iter->second.value;
        }  
    }
 
}

/**
 * @brief 
 * 
 * @param file 
 * @param n_step 
 * @param min_val 
 * @param max_val 
 */
void individual::output(std::ofstream& outfile, double n_step, double min_val, double max_val){
    double keep = fit_val;
    double env_value;
    
	for (int i = 0; i<n_step; i++){
        env_value = min_val + (i/n_step) * (max_val - min_val);
        fit_value(env_value);
		outfile << fit_val << "," ;
	}
    outfile << "\n";
    fit_val = keep;
};


population::population(int pop_size){
    ind = vector<individual>(pop_size);
}

population::~population(){

}

/**
 * @brief Goes to the next generation by calculate the fitness of each individual in the environment,
 * then sample in a poisson distribution for each individual to have their offspring number, finally add mutation
 * to the new individuals.
 * 
 *  //!// FOR NOW the fitness corresponds to the lambda parameter of the poisson for the number of offspring
 * 
 * @param env_value The current environmental value
 * @param param The parameters of the individual including the mutation rate
 */
void population::new_generation(double env_value, para_ind param){

    // RNG
    std::random_device rd;
    std::mt19937 gen(rd());

    // Calculate the fitness of each individual and add the sampled number of offspring in the 
    // new population.

    int K = 1000; //!// to add as parameter the limit capacity of the environment
    vector<individual> new_pop; // the vector that will contain the new generation

    vector<double> pop_fit; // the vector of fitness of each individual
    double tot_fit = 0; // The total fitness of the mutation considering the 
    int n_offspring = 0;

    for(int i = 0; i < ind.size(); i++){
        ind[i].fit_value(env_value);
        pop_fit.push_back(ind[i].fit_val);
        /* CHANGE POPULATION SIZE, NOT DONE
        //!// n_offspring = std::poisson_distribution<>(ind[i].fit_val)(gen); 

        for (int j = 0; j<n_offspring; j++){ //!// MAYBE CAN REMOVE THIS FORLOOP
            new_pop.push_back(ind[i]);
        }*/
    }
    std::discrete_distribution<>pop_distr(pop_fit.begin(), pop_fit.end());
    int random_integer;
     for(int i = 0; i < ind.size(); i++){
        random_integer = pop_distr(gen);
        new_pop.push_back(ind[random_integer]); //!// in that case do not have to use a pushback
    }

    /* SAME ONLY IF POPULATION SIZE CHANGE
    if (new_pop.size() > K){
        int size = new_pop.size();

        for (int i = 0; i < (size - K); i++){
            std::uniform_int_distribution<int> uni(0, new_pop.size()); 
            auto random_integer = uni(gen); 

            ind.erase(ind.begin() + random_integer);
            cout << random_integer << endl; //!//
        }
    }
    */

    ind = new_pop;
    int n_mut;
   
    for (int i = 0; i<ind.size(); i++){
        n_mut =  param.distr_mut_rat(gen);

        for (int m = 0; m < n_mut; m++){
            ind[i].chr.add_mutation(param);
        }
    }
}

void population::output_one_step(std::ofstream& outfile){
    for (int i = 0; i<1000; i++){
        outfile << i << "\t";
	}
    outfile << "\n";
	for (int i = 0; i<ind.size(); i++){
        ind[i].output(outfile);
	}
};

/**
 * @brief careful change the fit_val of the individual
 * 
 * @param outfile 
 */
void population::output_mean(std::ofstream& outfile){
    double n_step = 1000;
    double min_val = 0;
    double max_val = 1;
    double env_value;
    double sum = 0;


    vector<double> pop(ind.size(), 0);
    for (int i = 0; i<n_step; i++){
	    for (int j = 0; j<ind.size(); j++){
            env_value = min_val + (i/n_step) * (max_val - min_val);
            ind[j].fit_value(env_value);
		    pop[j] = ind[j].fit_val;
	    }
        for (int j = 0; j<ind.size(); j++){ 
            sum += pop[j];
        }
        outfile << sum / ind.size() << ",";
        sum = 0;
	}
    outfile << "\n";
};
