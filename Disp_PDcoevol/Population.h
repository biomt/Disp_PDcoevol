#pragma once
#include"Individuals.h"
#include <vector>

class Population {
public:
	Population();
	~Population();

	int N, Noff, Noff2; // counter of overal number of offsrping
	int Nf, Nm, bd_Nf,bd_Nm; // count of adults per sex and before dispersal numbers
	int Foff, Moff, Foff2, Moff2; // counter in order of how offspring is passed from through form offspring production to adulthood

	double mean_display,mean_pref, sd_display, sd_pref;
	
	
	std::vector<Individuals> females, males, femaleOffspring, maleOffspring, tmp_females, tmp_males, tmp2_females, tmp2_males;
	
	void outpop(int repl, int gen, int x, int y, double K, double S, std::ofstream *out);

};