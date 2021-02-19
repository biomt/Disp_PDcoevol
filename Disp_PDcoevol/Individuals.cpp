#include"Individuals.h"

// random number generator for the class individual, as the script structure is hierachical and can_t use the one in main
random_device rdInd;
mt19937 rdgenInd(rdInd());

Individuals::Individuals(bool sx, int coordx, int coordy) {

	sex = sx;
	x = coordx;
	y = coordy;

}

Individuals::~Individuals() {

}

void Individuals::setgenes(int L, std::normal_distribution<>NormPref, std::normal_distribution<>NormDisplay, std::normal_distribution<>NormEmig, std::normal_distribution<>NormEmigM) // this function is automatically used by the ind. apply the values for the alleles
{
	// one can use a more quant gen approach with assuming that the trait is additive. normal distribution of alles ends up
	// 

	for (int i = 0; i < L; i++)
	{
		traits.display1.push_back(NormDisplay(rdgenInd)); // fill chromosome up with alleles for loci drwan from a random gen
		traits.display2.push_back(NormDisplay(rdgenInd));

		traits.pref1.push_back(NormPref(rdgenInd));
		traits.pref2.push_back(NormPref(rdgenInd));

		traits.emig1.push_back(NormEmig(rdgenInd));
		traits.emig2.push_back(NormEmig(rdgenInd));

		traits.emigM1.push_back(NormDisplay(rdgenInd)); // error ??? why does it draw values different from emig1
		traits.emigM2.push_back(NormDisplay(rdgenInd));


		traits.g_display += traits.display1[i] + traits.display2[i]; // create additive trait value with equal effects of both chromosomes
		traits.g_pref += traits.pref1[i] + traits.pref2[i];
		traits.g_emig += traits.emig1[i] + traits.emig2[i];
		traits.g_emigM += traits.emigM1[i] + traits.emigM2[i];
		
	}
	//cout << "gemigM = " << traits.g_emigM << endl;
	

	if (traits.g_display > 0.0)
	{
	traits.p_display = traits.g_display; // phenotype equals the genotype when trait is above 0
	}
	
	else
	{
		traits.p_display = 0.0; // phenotype is bound by 0
	}


	traits.p_pref = traits.g_pref;

	if (traits.g_emig < 0.0)
	{
		traits.p_emig = 0.0; // emigration phenotype bound by 0
	}
	else 
	{
		if (traits.g_emig > 1.0)
		{
			traits.p_emig = 1.0; // emigration phenotype bound by 1
		}
		else {
			traits.p_emig = traits.g_emig;
		}
	}
	
	// male emig 
	if (traits.g_emigM < 0.0)
	{
		traits.p_emigM = 0.0; // male emigration phenotype bound by 0
	}
	else
	{
		if (traits.g_emigM > 1.0)
		{
			traits.p_emigM = 1.0; //male  emigration phenotype bound by 1
		}
		else {
			traits.p_emigM = traits.g_emigM;
		}
	}

}

// mutation effect missing emigration trait , individual output missing emig trait
void Individuals::mutation_effect(int L,int Nloci,normal_distribution<>effect)
{
	double muteffect;
	muteffect = effect(rdgenInd);

	

	if (L < Nloci * 2)
	{
		if (L < Nloci)
		{
			traits.pref1[L] += muteffect;
		}
		else {
			traits.pref2[L-Nloci] += muteffect; // 
		}

		traits.g_pref += muteffect;
		traits.p_pref = traits.g_pref;
	}
	else
	{
		if (L < Nloci * 4) 
		{

			if (L < Nloci * 3)
			{
				traits.display1[L - 2 * Nloci] += muteffect;
			}
			else {
				traits.display2[L - 3 * Nloci] += muteffect;	
			}

				traits.g_display += muteffect;
		
			if (traits.g_display > 0)
			{
				traits.p_display = traits.g_display; // phenotype equals the genotype when trait is above 0
			}

			else
			{
			traits.p_display = 0; // phenotype is bound by 0
			}


		}
		else
		{
			if (L < Nloci * 6) 
			{ //

				if (L < Nloci * 5)
				{
				traits.emig1[L - 4 * Nloci] += muteffect;
				}
				else {
				traits.emig2[L - 5 * Nloci] += muteffect;
				}

				traits.g_emig += muteffect;


				if (traits.g_emig < 0.0) traits.p_emig = 0.0;
				else
				{
					if (traits.g_emig > 1.0)traits.p_emig = 1.0;
					else traits.p_emig = traits.g_emig;

				}

			}
			else
			{
				if (L < Nloci * 7)
				{
					traits.emigM1[L - 6 * Nloci] += muteffect;
				}
				else {
					traits.emigM2[L - 7 * Nloci] += muteffect;
				}

				traits.g_emigM += muteffect;


				if (traits.g_emigM < 0.0) traits.p_emigM = 0.0;
				else
				{
					if (traits.g_emigM > 1.0)traits.p_emigM = 1.0;
					else traits.p_emigM = traits.g_emigM;

				}

			}

		}

		
	}
}

void Individuals::outind(int repl, int gen, std::ofstream *out)
{
	
	*out << repl << "\t" << gen << "\t" << x << "\t" << y << "\t" << sex << "\t" << traits.g_pref << "\t" << traits.p_pref << "\t" << traits.g_display << "\t" << traits.p_display << "\t" << traits.g_emig << "\t" << traits.p_emig << "\t" << traits.g_emigM << "\t" << traits.p_emigM << endl;

}

void Individuals::outind_slim(int repl, int gen, std::ofstream *out)
{

	*out << repl << "\t" << gen << "\t" << x << "\t" << y << "\t"  << traits.g_pref << "\t" << traits.g_display << "\t"<< traits.g_emig <<"\t" << traits.g_emigM << "\t" << sex << endl;

}

void Individuals::outind_dispersal(int repl, int gen, std::ofstream *out)
{

	*out << repl << "\t" << gen << "\t" << x << "\t" << y <<"\t"<< newx<<"\t"<<newy<<"\t"<< sex << "\t" <<dispersal_type<<"\t"<<traits.g_pref << "\t" << traits.g_display <<"\t" <<traits.g_emig << "\t" << traits.p_emig <<endl;

}





