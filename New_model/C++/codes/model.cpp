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
 * TODO: add a strength of selection rather than just fitness corresponding to probability to have an offsrping
 * TODO: problem of range of niche and range of environmental variation
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

    double env_min = 0; // The minimal value of the environment
    double env_max = 1; // The maximal value of the environment
    double amplitude = 0.5; // the amplitude of the fluctuation of the environ
    double stochasticity = 0; // The stochasticity coefficient
    int period = 100; // The period of 1 cycle in number of generation

    int burn_gen = 500;  // The number of generation with no output (to reach equilibrium)
    int keep_gen = 100;  // The number of generation with ouput (after the burn)

    int pop_size = 10; // The population size
    double mut_rate = 2; // The mutation rate (lambda parameter of a poisson distribution)
    double mean_val_mut = -0.002; // The mean value of a mutation (effect on fitness)
    double var_val_mut = 0.1; // The variance of the mutation value (effect on fitness)

    // If trade_off is on, each mutation will lead to another mutation.
    // This mutation will have a different range and a value of -1/10 of the original value
    bool trade_off = 0; 

    // Initialise the parameters

    parameters param( // contains all the parameters see the parameter.h for more info
        amplitude, stochasticity, period, env_min, env_max, // environment
        mut_rate, mean_val_mut, var_val_mut, trade_off // individual
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
    env.initialise(param.env);
    env.output(file_env);
    population pop(pop_size);

    // Burn the first generations to reach an equilbrium
    for (int i = 0; i<burn_gen; i++){
        cout << "burnt: "  << i << " / " << burn_gen << endl;
        pop.new_generation(env.optimum[i], param.ind);
    }

    // Keep the equilibrium generations
    for (int i = 0; i<keep_gen; i++){
        cout << "kept: "  << i << " / " << keep_gen << endl;
        pop.new_generation(env.optimum[i+burn_gen], param.ind);

        //!// heavy compute, should be remove or optimize if possible,
        // since go through all mutations several times

        pop.output_mean(ind_outfile); 
    }

    std::map<double, mutation>::iterator iter;
    individual* ad_ind = &pop.ind[0];

// // // // // // // // // // //

    ind_outfile.close();
    env_outfile.close();

    cout << "\n \n DONE \n \n";
    return 0;
}





