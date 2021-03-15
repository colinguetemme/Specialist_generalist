
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

	std::random_device rd; //!//
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> runif(0,1);

	double position = runif(gen);

	if (mutations.find(position) == mutations.end()) { // No mutation at that position DO WE NEED IT //!//
		mutation mut;
		double boundary1 = runif(gen);
		double boundary2 = runif(gen);
		mut.min = min(boundary1, boundary2);
		mut.max = max(boundary1, boundary2);
		mut.value = para.distr_mut_value(gen);

		mutations[position] = mut;
		n_mut++;
	}
}
