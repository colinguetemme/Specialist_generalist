#pragma once

#include <stdio.h>
#include <stdlib.h>
//#include <io.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>
#include <bitset>

using namespace std; 

//Parameters
int t; //simulation time  

const int n = 4;
int ng = pow(2, n);//number of gametes

double u; //mutation rate
double s; //strength of selection
double r[n-1]; //recombination rates

int L; //period
double A; //amplitude of fluctuation
double d; //magnitude of stochasticity

//Variables
double opt; //optimum


double alpha[n]; //allelic values

vector<double> p; //gamete frequencies

struct rec_table {
	double jk_table[n][n];  
};
vector<rec_table> rec_tables; //gamete values


vector<std::string> gametes;
vector<double> genotypes;


// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());

std::uniform_real_distribution<> unif(0.0, 1.0);


//Functions
void initialisation(void);
void recombination(void);


