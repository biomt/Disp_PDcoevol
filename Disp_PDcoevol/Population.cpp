#include"Population.h"

Population::Population() {

	N = 0;
	Noff = 0;
	Nf = 0;
	Nm = 0;
	bd_Nm = 0;
	bd_Nf = 0;
	Foff = 0;
	Moff = 0;
	Noff2 = 0; // number of offsrping after selection
	Foff2 = 0;
	Moff2 = 0;

	mean_display = 0;
	mean_pref = 0;
	sd_display = 0;
	sd_pref = 0;
}

Population::~Population() {
	females.clear();
	males.clear();
	femaleOffspring.clear();
	maleOffspring.clear();
}

void Population::outpop(int repl, int gen, int x, int y, double K,double S, std::ofstream *out)// output population
{
	*out << repl << "\t" << gen << "\t" << x << "\t" << y << "\t" << K << "\t"<< S << "\t"<< Nf << "\t" << Nm << endl;
}