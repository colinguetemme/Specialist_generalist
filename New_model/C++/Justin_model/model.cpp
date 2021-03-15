/**
 * @file model.cpp
 * @author Colin Guetemme
 * @brief Main file of the project
 * @version 0.1
 * @date 2021-03-10
 * 
 * @copyright Copyright (c) 2021
 * 
 */

# include "model.h"

// RNG GENERATOR OF THE PROJECT 
std::random_device rd;
std::mt19937 gen(rd());

int main() {
    
    
    ////////////////
    // PARAMETERS //
    ////////////////

    double env_min = 0;
    double env_max = 1;
    double amplitude = 0.5;
    double stochasticity = 0;
    int period = 24;

    double mut_rate = 0.01;
    double mean_val_mut = -0.1;
    double var_val_mut = 0.1;

    int pop_size = 200;

    // Initialise the parameters

    parameters param( // contains all the parameters see the parameter.h for more info
        amplitude, stochasticity, period, env_min, env_max, // environment
        mut_rate, mean_val_mut, var_val_mut // individual
        );


    /////////////////
    // THE OUTPUTS //
    /////////////////

    ///////////////
    // THE MODEL //
    ///////////////

    environment env;
    env.initialise(param);

    population pop(pop_size);

    for (int i = 0; i<env.optimum.size(); i++){
        //cout << env.optimum[i] << endl;
        pop.new_generation(env.optimum[i], param.ind);
    }

    std::map<double, mutation>::iterator iter;
    individual* ad_ind = &pop.ind[0];

    /*
    for(iter = pop.ind[0].chr.mutations.begin(); iter != pop.ind[0].chr.mutations.end(); iter++){   
        cout << "mut " << iter->second.value << endl;
    }
    */ 

    for (int i = 0; i< pop.ind.size(); i++){
        cout << "individual" << i << " : " << pop.ind[i].fit_val << endl;
    }
    
    return 0;
}




