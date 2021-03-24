# include "parameters.h"

/**
 * @brief Construct a new parameters::parameters object
 * 
 * @param a The amplitude of the environment variation
 * @param d The stochasticity of the environment variation
 * @param L The period of the environment variation
 * @param min The minimum of the environment variation
 * @param max The maximum of the environment variation
 * 
 * @param mr The mutation rate of the population
 * @param mvm The mean value of a mutation
 * @param vvm The variance of the value of a mutation
 */
parameters::parameters(
    double a, double d, int L, double min, double max,
    double mr, double mvm, double vvm
    ){
    env.amplitude = a;
    env.stochastic = d;
    env.period = L;
    env.min = min;
    env.max = max;

    ind.mut_rate = mr;
    ind.mean_val_mut = mvm;
    ind.var_val_mut = vvm;
    ind.distr_mut_rat = std::poisson_distribution<>(ind.mut_rate);
    ind.distr_mut_value = std::normal_distribution<>(ind.mean_val_mut, ind.var_val_mut);
}

parameters::~parameters(){

}
