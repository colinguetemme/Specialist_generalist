/**
 * @file model.cpp
 * @author Colin Guetemme
 * 
 * @brief This model was described by Justin Travis, the basic idea is to have a population where each 
 * individual as a continuous chromosome. 
 * The generation does not overlap and at each generation we sample with replacement in the previous
 * generation with a probability based on the fitness of each individual.
 * The fitness of each individual depends of an environment which changes through time 
 * (following a sinus function).
 * From one generation to another, individual can get mutations, each mutation having a range of activation
 * (e.g. from pH 2 to 3 if using pH) and a value, which represents the effect of fitness of this mutation.
 * To know the fitness of an individual for a given parameter value, we do the some of all the mutations for
 * which the range include the environmental parameter.
 * 
 * TODO: Every x generations output the all phenotype
 * TODO: red noise (positively autocorrelated with parameter l)
 * TODO: mutpop, keep track of each mutation in the simulation, maybe just the end
 * TODO: add a non constant population model (just a K model)
 * TODO: add expression level (each value of the mutation multiply by expression level)
 * 
 * The parts of the code that need particular attention (bugs, optimization, ...)
 * are commented with //!//
 * 
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

    // TO DECLARE WHICH TYPE OF ENVIRONMENT TO USE
    // options: "sinus", "custom"
    string env_type = "sinus";

    double env_min = 0; // The minimal value of the environment
    double env_max = 1; // The maximal value of the environment
    double amplitude = 0.95; // the amplitude of the fluctuation of the environ
    double stochasticity = 0; // The stochasticity coefficient //!// CHECK FOR CUSTOM
    int period = 100; // The period of 1 cycle in number of generation

    vector<double> custom_states {0.1, 0.9}; // 

    int burn_gen = 0;  // The number of generation with no output (to reach equilibrium)
    int keep_gen = 1000;  // The number of generation with ouput (after the burn)
    int step_keep = 5; // Output one generation out of step_keep generation (e.g. 1 every 5)
    int tmax = burn_gen + keep_gen;

    int pop_size = 100; // The population size
    double mut_rate = 0.01; // The mutation rate (lambda parameter of a poisson distribution)
    double prop_del = 0.5; // The proportion of deleterious mutation (compare to beneficial)
    double mean_val_mut = 0.05; // The mean value of a mutation (effect on fitness)
    double var_val_mut = 0.01; // The variance of the mutation value (effect on fitness)
    //!// Deleterious should be 1000 more frequent than beneficial
    bool uni = 1; //!// in test, if = 1 all the value of parameter are equally likely to be in a mutation
    double min_fit = 0.05; //!// not added, minimal value of the fitness to 

    // If trade_off is on, each mutation will lead to another mutation.
    // This mutation will have a different range and a value of -1/2 of the original value
    bool trade_off = 0; 
    bool trade_off_del = 0; // also trade-off for deleterious
    double trade_off_strength = 0.3; // The strength of the trade off (proportion of the original value)

    // Initialise the parameters

    parameters param( // contains all the parameters see the parameter.h for more info
        amplitude, stochasticity, period, env_min, env_max, tmax, // environment
        mut_rate, mean_val_mut, var_val_mut, prop_del, trade_off, trade_off_del, trade_off_strength // individual
        );

    /////////////////
    // THE OUTPUTS //
    /////////////////
 
    // Creates the dir for that simulation 
    string dirname = "../results/data/";

    // Creates the file to output the environment
    string file_env = dirname + "env.csv";
    std::ofstream env_outfile (file_env);

    // Creates the file to output the population results
    string file_ind = dirname + "ind.csv";
    std::ofstream ind_outfile (file_ind);

    ///////////////
    // THE MODEL //
    ///////////////

    environment env;

    if (env_type == "sinus"){
        env.initialise_sinus(param.env);
    } else if (env_type == "custom"){
        env.initialise_custom(param.env, custom_states);
    }
     
    env.output(file_env);
    population pop(pop_size);

    // Burn the first generations to reach an equilbrium
    for (int i = 0; i<burn_gen; i++){
        cout << "burnt: "  << i << " / " << burn_gen << endl;
        //pop.new_generation(env.optimum[i], param.ind);
        pop.new_generation_K(env.optimum[i], param.ind, 10000); //!// not working change the population size
    }

    // Keep the equilibrium generations
    for (int i = 0; i<keep_gen; i++){
        
        //pop.new_generation(env.optimum[i], param.ind);
        pop.new_generation_K(env.optimum[i], param.ind, 10000); //!// not working change the population size

        //!// heavy compute, should be minimize the number of iteration here,
        // since go through all mutations several times
        if (i % step_keep == 0){
            cout << "kept: "  << i << " / " << keep_gen << endl;
            pop.output_mean(ind_outfile); 
        }
    }

    std::map<double, mutation>::iterator iter;
    individual* ad_ind = &pop.ind[0];

// // // // // // // // // // //

    ind_outfile.close();
    env_outfile.close();

    cout << "\n DONE \n \n";
    return 0;
}


