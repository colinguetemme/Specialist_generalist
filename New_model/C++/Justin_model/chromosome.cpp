#include  "chromosome.h"


Chromosome::Chromosome() {
	nMut = 0;
	Nho = 0;
}


Chromosome::~Chromosome() {
}

///////////////
/// METHODS ///
///////////////


void Chromosome::deleteChromo(void) {
	if (!mutations.empty()) mutations.clear();
}



double Chromosome::add_mutation(int homolog, double position, double selection, double dominance) {
	double viability;
	mutation mut;
	mut.homol = homolog;
	mut.s = selection; 
	mut.h = dominance;

	if (mutations.find(position) == mutations.end()) { // No mutation at that position
		nMut++;
		mutations[position] = mut;
		viability = (1.0 - mut.h * mut.s);
	}
	else viability = 1.0;

	return viability;
}
