
#pragma once

///////////////
// LIBRARIES //
///////////////

#include <map>
#include <set>

#include <stdio.h>
#include <fstream>
#include <iostream>

//#include "Parameters.h"

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
//typedef std::map<double, int, std::less<double>> MapPos;

//map containing continuous allele coding for a trait: map<position, allele>
//typedef std::map<double, double, std::less<double>> MapTrait;


/**
 * @brief hou
 * 
 */
class Chromosome {
public:

    // Constructor / Destructor
	Chromosome();
	~Chromosome();

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
     * @param homolog if the mutation is homologue (0 on 1, 1 on 1, 2 on 2)
     * @param position position in the continuous chromosome
     * @param selection not used
     * @param dominance not used
     * 
     * @return The viability of the Chromosome
     * 
     */
	double add_mutation( 
		int, // homologue
		double, // position
		double, // s (selection coefficient)
		double // h (dominance coefficient)
		); 

    /**
    * @brief Delete all the mutations on a chromosome
    * 
    */
	void deleteChromo();

};