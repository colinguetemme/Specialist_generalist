
#pragma once

///////////////
// LIBRARIES //
///////////////

#include <map>
#include <set>
#include <random>

#include <stdio.h>
#include <fstream>
#include <iostream>

#include "parameters.h"

using namespace std;

////////////////
// DEFINITION //
////////////////

/**
 * @brief Contains the range of effect of the mutation and the effect on fitness of this mutation
 * 
 */
struct mutation {

	double min; // Minimum of the range of effect of the mutation
	double max; // Maximum of the range of effect of the mutation

	double value; // Effect on the fitness
};

// Map containing all the mutations
typedef std::map<double, mutation, std::less<double>> MapMuts; 

// Map containing all the mutation position
typedef std::map<double, int, std::less<double>> MapPos;

//map containing continuous allele coding for a trait: map<position, allele>
//typedef std::map<double, double, std::less<double>> MapTrait;


/**
 * @brief Class describing the continuous chromosome and 
 * 
 */
class chromosome {
public:

    // Constructor / Destructor
	chromosome();
	~chromosome();

    ///////////////////
	// THE VARIABLES //
    ///////////////////

	int n_mut; // The number of mutations on that chromosome
	MapMuts mutations; // The map with the position of each mutation on the chromosome (see MapMuts in chromosme.h)

    /////////////////
	// The METHODS //
    /////////////////
    
    /** @brief Add a mutation to the continuous chromosome, 
     * the mutation parameters will be randomly sampled when the function is called
     *
     * @param param the structure containing all the parameters of the individual
     */
	void add_mutation(para_ind); 
    void add_mutation2(para_ind); 

    /**
    * @brief Delete all the mutations on a chromosome
    * 
    */
	void delete_chromo();

};
