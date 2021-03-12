# include "parameters.h"

parameters::parameters(double a, double d, int L){
    env.amplitude = a;
    env.stochastic = d;
    env.period = L;
}

parameters::~parameters(){

}

para_env::para_env(){

};
para_env::~para_env(){
    
};