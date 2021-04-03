# include "parameters.h"


/**
 * @brief Construct a new parameters::parameters object
 * 
 * @param a The amplitude of the environment variation
 * @param d The stochasticity of the environment variation
 * @param L The period of the environment variation
 * @param min The minimum of the environment variation //!// not useful
 * @param max The maximum of the environment variation //!// not useful
 * @param tmax The maximum duration of the simulation
 * 
 * @param mr The mutation rate of the population
 * @param mvm The mean value of a mutation
 * @param vvm The variance of the value of a mutation
 * @param p_d The proportion of deleterious in the population
 * @param t_o If 1, each beneficial mutation as a deleterious trade-off
 * @param t_o_d If 1, deleterious as also a beneficial trade-off
 * @param t_o_s The strength of the trade-off (=the proportion of the value of the original mutation)
 */
parameters::parameters(
    double a, double d, int L, double min, double max, int tmax,
    double mr, double mvm, double vvm, double p_d, bool t_o, bool t_o_d, double t_o_s
    ){
    env.amplitude = a;
    env.stochastic = d;
    env.period = L;
    env.min = min;
    env.max = max;
    env.tmax = tmax;

    ind.mut_rate = mr;
    ind.mean_val_mut = mvm;
    ind.var_val_mut = vvm;
    ind.trade_off = t_o;
    ind.prop_del = p_d;
    ind.trade_off_del = t_o_d;
    ind.trade_off_strength = t_o_s; 
    ind.distr_mut_rat = std::poisson_distribution<>(ind.mut_rate);
    ind.distr_mut_value = std::normal_distribution<>(ind.mean_val_mut, ind.var_val_mut);
     
}

parameters::~parameters(){

}
