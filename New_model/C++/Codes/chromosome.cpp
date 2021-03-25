
#include  "chromosome.h"


chromosome::chromosome() {
}

chromosome::~chromosome() {
}

///////////////
/// METHODS ///
///////////////


void chromosome::delete_chromo(void) {
	if (!mutations.empty()) mutations.clear();
}


/**
 * @brief Add randomly a mutation on the chromosome defined as a map
 * 
 * @param para the parameters of the individuals
 */
void chromosome::add_mutation(para_ind para) {

	// Create the random device 
	std::random_device rd; //!// 
	std::mt19937 gen(rd());

	// The require distribution for the mutation's range of effect
	std::uniform_real_distribution<> runif(-3,4);

	double position = runif(gen);

	if (mutations.find(position) == mutations.end()) { // No mutation at that position DO WE NEED IT //!//
		mutation mut;

		// Start by sampling two values from the distribution and then look at which one is the min and max
		double boundary1 = runif(gen);
		double boundary2 = runif(gen);
		mut.min = min(boundary1, boundary2);
		mut.max = max(boundary1, boundary2);

		// The value of the mutation is sampled from the corresponding normal distribution
		mut.value = para.distr_mut_value(gen);

		// Add the mutation on the map of mutation of that chromosome
		mutations[position] = mut;


		// Add a second mutation to generate the trade-off of this mutation
		// is not counted in the number of mutation (because only one)
		if (para.trade_off == 1){
			double position_t_o = runif(gen); // _t_o stands for trade-off
			if (mutations.find(position_t_o) == mutations.end()) {
				mutation mut_t_o;

				double boundary1 = runif(gen);
				double boundary2 = runif(gen);
				mut_t_o.min = min(boundary1, boundary2);
				mut_t_o.max = max(boundary1, boundary2);
				mut_t_o.value = -mut.value/10; //!// should be a parameter and maybe different if positive and negative
				mutations[position_t_o] = mut_t_o;
			}

			n_mut++; // needed ? //!//
		}
	}
}
