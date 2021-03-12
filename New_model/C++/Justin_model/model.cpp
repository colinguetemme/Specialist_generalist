/**
 * @file model.cpp
 * @author Colin Guetemme (you@domain.com)
 * @brief Main file of the project
 * @version 0.1
 * @date 2021-03-10
 * 
 * @copyright Copyright (c) 2021
 * 
 */

using namespace std;

# include "parameters.h"
# include "model.h"

#include <iostream>

int main() {
    double amplitude = 0.5;
    double stochasticity = 0.5;
    int period = 24;
    parameters param(amplitude, stochasticity, period);
    cout << param.env.amplitude << endl;
    return 0;
}

