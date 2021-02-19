#pragma once
#include<vector>
using namespace std;


class trait{

public:
	
	trait();
	~trait();

	vector<double> pref1, pref2; // homologue chromosomes with loci of the trait
	vector<double> display1, display2;
	vector<double> emig1, emig2;
	vector<double> emigM1, emigM2;

	// variables

	double viability;

	// genotypic value
	
	double g_pref, g_display,g_emig,g_emigM;

	// phenotypic value

	double p_pref, p_display,p_emig,p_emigM;

	

};

// possible exercise, recombination by having loci in one vector Greta style