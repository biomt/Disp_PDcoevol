#include "SpatialRunaway.h"

// DEFINE DISTRIBUTIONS












#if CLUSTER

int main(int argc,char* argv[])
{
	//Get the curretn directroy,  Always run in every model. Saves path of the current directory and put input in there
	char*buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; // current directory path
	dirout = dir + "Outputs/"; 

	simNr = std::atoi(argv[1]); // the column for which parameters to change when running in the cluster
	rep = std::atoi(argv[2]);
	gen = std::atoi(argv[3]);

	sample_all_mates=std::atoi(argv[4]);
	ind_sample=std::atoi(argv[5]);
	quantileS=std::atoi(argv[6]);
	subsample_mates=std::atoi(argv[7]);

	disprob=std::atof(argv[8]);
	
	w_m=std::atoi(argv[9]);
	f_costs=std::atoi(argv[10]);
	w_f=std::atoi(argv[11]);

	Kmin=std::atoi(argv[12]);
	Kmax=std::atoi(argv[13]);

	Smin=std::atoi(argv[14]);
	Smax=std::atoi(argv[15]);

	out_pop_interval=std::atoi(argv[16]);
	out_ind_interval=std::atoi(argv[17]);
	out2_limit=std::atoi(argv[18]);

	out_pop_interval2=std::atoi(argv[19]);
	out_ind_interval2=std::atoi(argv[20]);

//addition for dispersal evolution
	disp_evol= std::atoi(argv[21]);
	disp_sex= std::atoi(argv[22]);
	disp_cost= std::atof(argv[23]);

	meanDisplay= std::atof(argv[24]);
	meanPref= std::atof(argv[25]);
	meanEmig= std::atof(argv[26]);
	meanEmigM = std::atof(argv[27]);

	stdDisplay= std::atof(argv[28]);
	stdPref= std::atof(argv[29]);
	stdEmig= std::atof(argv[30]);
	stdEmigM = std::atof(argv[31]);

	
// track time of the line

	extime = clock();

	runmodel();
	//cin.get();// stop befor finishing and keep console open
	return 0;

}


#else
int main(void)
{
	//Get the curretn directroy,  Always run in every model. Saves path of the current directory and put input in there
	char*buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; // current directory path
	dirout = "C:\\Users\\r01mt18\\Documents\\Models\\dispersal_PDcoevol\\Outputs\\";
	// track time of the line

	extime = clock();

	runmodel();
	cin.get();// stop befor finishing and keep console open
	return 0;

}

#endif



void runmodel(void){

	out_pop_header();// creates the files
	out_ind_header();
	if(dispersal_out) out_inddispersal_header();
	out_param_header();
	outparam(simNr,&param);

	



	for (r = 0; r < rep; r++) // the replicate loop for your model
	{
		cout << "replicate " << r << endl;
		
		start();				// initilise each 

		//cout << "start OK" << endl;
		
		for (g = 0; g < gen; g++) // generation loop !!!!!!!!!!!!!!
		{
			cout << "generation " << g << endl;
			
			sim_out(); //  output 

			reproduction2();
			//cout << "reproduction OK" << endl;
			
			if (d_post) { // post-survival dispersal

				survival();
				//cout << "survival OK" << endl;

				if (dispersal_cntrl)
				{
					dispersal_control();
				}
				else
				{
					dispersal();
				}
			}
			else { // pre-survival dispersal

				dispersal_ds();
				survival_ds();

			}
			
			
			//cout << "dispersal OK" << endl;
		}

		for (int i = 0; i < xmax; i++) {
			for (int j = 0; j < ymax; j++) {
				popgrid[i][j].females.clear();
				popgrid[i][j].males.clear();
				popgrid[i][j].N = 0;
				popgrid[i][j].Nm = 0;
				popgrid[i][j].Nf = 0;
			}
		}
	}

	pops.close();
	inds.close();
	if(dispersal_out) inds_dispersal.close();
	param.close();
	

}// masterfunction of your model, calls all your functions

void start(void) {

	normal_distribution<> NormPref(meanPref / (2.0*(double)Nloci), stdPref / sqrt(2.0*(double)Nloci));
	normal_distribution<> NormDisplay(meanDisplay / (2.0*(double)Nloci), stdDisplay / sqrt(2.0*(double)Nloci));
	normal_distribution<> NormEmig(meanEmig / (2.0*(double)Nloci), stdEmig / sqrt(2.0*(double)Nloci));
	normal_distribution<> NormEmigM(meanEmigM / (2.0*(double)Nloci), stdEmigM / sqrt(2.0*(double)Nloci));

	uniform_int_distribution<> Kdistr(Kmin, Kmax); // distr for carrying capacity
	uniform_int_distribution<> Sdistr(Smin, Smax); // percent of males subsampled in a population

	for (int i = 0; i < xmax; i++)
	{
		for (int j = 0; j < ymax; j++)
		{
			grid[i][j].K = Kdistr(rdgen); // samples a value from my uniform distribution and applies it to the cell
			popgrid[i][j].N = (int)grid[i][j].K; // set the carrying capacity for this population

			if (sample_all_mates){}

			else
			{
				grid[i][j].S = Sdistr(rdgen); // set

				
				//while ((int)grid[i][j].S % 10 != 0) this while statment would be necesarry if I would want percents 1-100 
				{
					//grid[i][j].S = Sdistr(rdgen);
				}

				//cout <<"S is "<< grid[i][j].S << endl;
			}


			for (int l = 0; l < (int)(popgrid[i][j].N / 2); l++) // intilise population with half females and half males
			{
				popgrid[i][j].females.push_back(Individuals(1, i, j)); // in the cell, initialise the females of the pop
				popgrid[i][j].females[l].setgenes(Nloci, NormPref, NormDisplay,NormEmig,NormEmigM); // give the individuals their properties, genome..
				popgrid[i][j].Nf++;
			}

			for (int l = 0; l < (int)(popgrid[i][j].N) / 2; l++)
			{
				popgrid[i][j].males.push_back(Individuals(0, i, j)); // do the same for the males.
				popgrid[i][j].males[l].setgenes(Nloci, NormPref, NormDisplay, NormEmig,NormEmigM);
				popgrid[i][j].Nm++;
			}

		}


	}
}


void reproduction2(void) 
{
	bernoulli_distribution extdistr(hab_turn); // distribution of population exticntion, probability distribution of being true.
	uniform_real_distribution<> Prob(0.0, 1.0);
	poisson_distribution<>offdistr(offspring);
	bernoulli_distribution sexdistr(0.5); // prob of true (1) , could change the sex ratio with that

	// define local variable for the funtion
	vector<Individuals>::iterator iter; // member of the vector class, enables you to act on the vector
	vector<Individuals>::iterator iter2; // to loop through all males to calculate the sum of exponential pref to display relation
	vector<Individuals>::iterator iter3; // to loop through each male to apply cumultative distribution

	int IDmale;
	double sumNmales = 0.0; // variable for the sum of all males realtion of their dispaly to the females pref
	double probmale = 0.0; // probability of each male
	double cumprob = 0.0; // cumultative prob,
	double choice;

	int a_males; // counter of males
	int noffspring;
	bool sex; // sex of the offspring

	Individuals *ind; // is a pointer to now an undfeined object of class individual

	vector<double>cumdistr; // is the cumulative distribution of males prob of being chosen
	vector<int>idmates;	// identity of the mates, position of the mates in the male vector
	vector<int>a_mates; // males to sample from if not all males

	for (int i = 0; i < xmax; i++) // looping through the cell grid again
	{
		for (int j = 0; j < ymax; j++)
		{
			//deciding if this grid cell is going extinct

			if (extdistr(rdgen))
			{


				if (!popgrid[i][j].females.empty())
				{
					popgrid[i][j].females.clear();
					popgrid[i][j].Nf = 0;
				}


				if (!popgrid[i][j].males.empty())
				{
					popgrid[i][j].males.clear();
					popgrid[i][j].Nm = 0;
				}

				popgrid[i][j].N = 0;


			}

			else {}





			if (popgrid[i][j].Nm > 0 && popgrid[i][j].Nf > 0) // if there are males and females to reproduce
			{
				//cout << " R:n " << popgrid[i][j].N << endl;
				//cout << " R:nfemales " << popgrid[i][j].Nf << endl;
				//cout << "R: nmales " << popgrid[i][j].Nm << endl;

				for (iter = popgrid[i][j].females.begin(); iter != popgrid[i][j].females.end(); iter++)// use iterator to loop through female vector. For each female describe how she gets her mate:
				{
					sumNmales = 0.0;// initilalize at zero for every female
					probmale = 0.0;
					cumprob = 0.0;

					//cout << "female pref value: " << iter->traits.p_pref << endl;

					if (sample_all_mates || (!sample_all_mates && ind_sample && popgrid[i][j].Nm <= subsample_mates)) // dont have to specify further, cause if its boolean its automaitcally true for if statments. If all males should be sampled. 
					{
						// 1. get cumultative probabilty distribution for females to sample from, fitst loop through all males and calculate sum of probabiltiy
						for (iter2 = popgrid[i][j].males.begin(); iter2 != popgrid[i][j].males.end(); iter2++)
						{
							sumNmales += exp(iter2->traits.p_display*iter->traits.p_pref); // calculate sum, in the denominator of the equation for calcualting prob for each male!!!!!!!!!!!!!!!!!!!!!!!!!!!!

							//double pippo = exp(iter2->traits.p_display*iter->traits.p_pref);
							//cout <<"pippos should add to sumNmales"<< pippo << endl;

						}

						//cout  << "sumNmales: "<< sumNmales << endl;

						for (iter3 = popgrid[i][j].males.begin(); iter3 != popgrid[i][j].males.end(); iter3++) // calculate prob for each male and create cumultative distribution
						{
							probmale = (exp(iter3->traits.p_display*iter->traits.p_pref) / sumNmales);
							cumprob += probmale;
							cumdistr.push_back(cumprob);


							//cout <<" .male probability of being chosen : " << probmale << " cumprob: " << cumprob << " male display trait: "<< iter3->traits.p_display<< endl;

						}

						//2. let female choose : 

						//random number between 0-1
						choice = Prob(rdgen);

						//cout << "choice" << choice << endl;

						IDmale = 0;  // 

						//int pippo = (int)cumdistr.size();
						//out << pippo << endl;

						for (int z = 0; z < (int)cumdistr.size(); ++z) // choosing the female by iterating of the vector that contains the added probabilities, and therefore acts as a cumultative discribution). bug in the choice!!!
						{
							if (choice < cumdistr[z])
							{
								IDmale = z;
								break;
							}

						}

						//cout << "male ID " << IDmale << endl;

						cumdistr.clear();



					}
					
					else // if not sample all males
					{
						//cout << "subsample" << endl;

						if (ind_sample) // if a specific number of males is to be sampled
						{


							//sampling without  replacement

							a_males = popgrid[i][j].Nm;//the number of so called available males to take a sample from

							for (int z = 0; z < a_males; ++z) // gives males a ID number according to their place in the vecor a_males
							{
								a_mates.push_back(z);

								//cout << a_males<<  endl;
							}

							//cout << subsample_mates << endl;

							for (int z = 0; z < subsample_mates; ++z)
							{
								uniform_int_distribution<>samplemates(0, a_males - 1);

							

								IDmale = samplemates(rdgen);

								//cout << IDmale << "position of sampled male in availavle male vector "<< endl;

								sumNmales += exp(popgrid[i][j].males[IDmale].traits.p_display*iter->traits.p_pref);// sum of male prob


								idmates.push_back(IDmale);
								a_males--;
								a_mates.erase(a_mates.begin() + IDmale);

								
								//cout <<"Nm is" << popgrid[i][j].Nm<<" and subsample number is" <<sample <<" and subsample INT is" << (int)sample << endl;

							}


							if (!a_mates.empty()) a_mates.clear();// if vector is not empty, clear it
						}

						else // this is if a percentage of the population is to be sampled
						{
							//cout << "percents" << endl;
							a_males = popgrid[i][j].Nm;//the number of so called available males to take a sample from
							double sample = 1.0;

							if (quantileS) 
							{
								sample = (static_cast<double>(a_males) / 10.0)*static_cast<double>(grid[i][j].S);
							}

							else
							{
								sample = (static_cast<double>(a_males) / 100.0)*static_cast<double>(grid[i][j].S);
							}

							

							if (sample < 1.0) sample = 1; // if number of males is below 1, then the number of subsampled males is always 1

							//int sample = ((int)a_males /(int)10)*(int)grid[i][j].S; if subsampling is grouped into steps of 10 percent
							//if (sample < (int)1) sample = a_males; // if number of males is below 10, then the number of subsampled males is always 1

							//cout << a_males<<"and"<< a_males/10<<"and"<< grid[i][j].S << endl;
							//cout <<"Nm is" << popgrid[i][j].Nm<<" and subsample number is" <<sample <<" and subsample INT is" << (int)sample << endl;

							for (int z = 0; z < a_males; ++z) // where is the vector a_mates vector?????????
							{
								a_mates.push_back(z);
							}

							for (int z = 0; z < (int)sample; ++z)
							{
								uniform_int_distribution<>samplemates(0, a_males - 1);

								IDmale = samplemates(rdgen);

								//cout << IDmale << "position of sampled male in availavle male vector " << endl;
								sumNmales += exp(popgrid[i][j].males[IDmale].traits.p_display*iter->traits.p_pref);// sum of male prob


								idmates.push_back(IDmale);
								a_males--;
								a_mates.erase(a_mates.begin() + IDmale);

							}


							if (!a_mates.empty()) a_mates.clear();// if vector is not empty, clear it

						}
						
						// SAMPLED MALES ARE IN IDMATES VECTOR- READY TO BE CHOSEN BY FEMALE
						for (int z = 0; z < (int)idmates.size(); ++z)
						{
							probmale = (exp(popgrid[i][j].males[idmates[z]].traits.p_display*iter->traits.p_pref) / sumNmales); // enter male value of pogrid i and j at position determined by idmates of z
							
							//cout << probmale << " probability of male" << endl;
							cumprob += probmale;
							cumdistr.push_back(cumprob);

						}

						//cout << cumprob << " cumlatiative prob " << endl;

						choice = Prob(rdgen);

						//cout << choice << " her choice" << endl;

						IDmale = 0;

						//cout << choice << endl;

						for (int z = 0; z < (int)idmates.size(); ++z) // idmates and cumdistr vector refer both to the same individuals, idmates holds the position int the male vector, while cumdistr his probability of being chosen
						{
							if (choice < cumdistr[z])
							{
								IDmale = idmates[z]; // now ID male carries the position of the chosen male in the male vector

								//cout << IDmale << " choice process" << endl;
								break;
							}
						}

						//cout << IDmale << " IDmale" << endl;

						if (!idmates.empty()) idmates.clear();//

						cumdistr.clear();

					}
					// female has a male identified in IDmale

					//cout << "mating OK" << endl;

					// producing offpsring with male

					noffspring = offdistr(rdgen);
					
					//noffspring =fmax(1, offdistr(rdgen)); // changes sinced versin october, now number of offspring is min 1 per female, and max drawn from poisson distribution.
					//cout << noffspring << endl;

					if (noffspring > 0)
					{
						for (int z = 0; z < noffspring; ++z)
						{
							sex = sexdistr(rdgen);

							if (sex)// female
							{
								ind = new Individuals(sex, i, j); // with the pointer i know the adress of the object without using a counter
								inheritance(ind, *iter, popgrid[i][j].males[IDmale]);

								popgrid[i][j].femaleOffspring.push_back(*ind);// object is accessed by the star
								popgrid[i][j].Foff++;

							}
							else {
								ind = new Individuals(sex, i, j);
								inheritance(ind, *iter, popgrid[i][j].males[IDmale]);

								popgrid[i][j].maleOffspring.push_back(*ind);
								popgrid[i][j].Moff++;

							}

							popgrid[i][j].Noff++;
							delete ind;


						}
					}

					else{}
					//cout << "offspring OK" << endl;



				}


				/////Population level in the loop. Put the mutation function
				
				if (popgrid[i][j].Noff > 0)
				{
					mutation(i, j); 
				}

				// aplly mutations on the pop level

				//cout << "mutation OK" << endl;

			}
			// killing of the adult population
			popgrid[i][j].females.clear();
			popgrid[i][j].males.clear();

			//cout << "R: N= " << popgrid[i][j].N << endl;
			//cout << "R: Noff= " << popgrid[i][j].Noff << endl;

			//cout << "R: Nf= " << popgrid[i][j].Nf << endl;
			//cout << "R: Foff= " << popgrid[i][j].Foff << endl;

			//cout << "R: Nm= " << popgrid[i][j].Nm << endl;
			//cout << "R: Moff= " << popgrid[i][j].Moff << endl;

			popgrid[i][j].N = 0;
			popgrid[i][j].Nf = 0;
			popgrid[i][j].Nm = 0;
		}

	}
}


void inheritance(Individuals *kid, Individuals mummy, Individuals daddy)
{ 
	// offspring chromosomes with loci pref1,display1 and emig1 are inherited from mummy
	bernoulli_distribution berndistr(0.5); // prob of true (1) general bernoulli
	bernoulli_distribution recomb(recprob); // probability of recombination
	

	int rdn, rdn2, rdn3, rdn4, rdn5, rdn6, rdn7,rdn8; // sample from which homolog to start for each pair 

	rdn = berndistr(rdgen); // starting with either of the mum's chromosomes for loci encoding pref
	rdn2 = berndistr(rdgen); // -||- dad pref
	rdn3 = berndistr(rdgen);//  -||- mum display	
	rdn4 = berndistr(rdgen); // -||- dad dispaly
	rdn5 = berndistr(rdgen); // mum emig
	rdn6 = berndistr(rdgen); // dad emig
	rdn7 = berndistr(rdgen); // mum emigM
	rdn8 = berndistr(rdgen); // dad emigM

	for (int z = 0; z < Nloci; ++z) // to minimise loops, we punch in all the chromosomes. Reasons to keep it seperate would be if a trait shall be fixed and cant evolve
	{
		// mother loci for preference, to set pref1 of offspring
		if (rdn) rdn -= recomb(rdgen); // decinding if this loci is a site of recombination, if rdn is 1 (which is maternal chromosome number 1) , then recombiantions happens by drawing from recomb by substratcting
		else rdn += recomb(rdgen);// if rdn is 0 (which is maternal chromosome2), then recombination happens with drawing from recomb and addition
		
		// 1 in this case does not represnt recombination,
		// instead 0 or 1 describes the position of the inherited loci at the parental chromosome 
		// and is the outcome of which chromosome was chosen to start with 
		// and if recombiantion happend subsequently

		if(rdn)//  (setting the gametes so to say)
		{ 
			kid->traits.pref1.push_back(mummy.traits.pref1[z]); // 
			kid->traits.g_pref += mummy.traits.pref1[z];

		}
		else{
			kid->traits.pref1.push_back(mummy.traits.pref2[z]);
			kid->traits.g_pref += mummy.traits.pref2[z];
		}

		 // dads loci for pref, to set pref2 of offpsring
		if (rdn2)rdn2 -= recomb(rdgen); // 
		else rdn2 += recomb(rdgen);

		if (rdn2)//
		{
			kid->traits.pref2.push_back(daddy.traits.pref1[z]);
			kid->traits.g_pref += daddy.traits.pref1[z];

		}
		else {
			kid->traits.pref2.push_back(daddy.traits.pref2[z]);
			kid->traits.g_pref += daddy.traits.pref2[z];
		}

		// mums loci for display, to set display1 of offspring
		if (rdn3)rdn3 -= recomb(rdgen); // 
		else rdn3 += recomb(rdgen);// 

		if (rdn3)// 
		{
			kid->traits.display1.push_back(mummy.traits.display1[z]); // 
			kid->traits.g_display += mummy.traits.display1[z]; // 

		}
		else {
			kid->traits.display1.push_back(mummy.traits.display2[z]);
			kid->traits.g_display += mummy.traits.display2[z];
		}

		// dads loci for display, to set display2 in offspring
		if (rdn4)rdn4 -= recomb(rdgen); // 
		else rdn4 += recomb(rdgen);// 

		if (rdn4)//
		{
			kid->traits.display2.push_back(daddy.traits.display1[z]); // 
			kid->traits.g_display += daddy.traits.display1[z];

		}
		else {
			kid->traits.display2.push_back(daddy.traits.display2[z]);
			kid->traits.g_display += daddy.traits.display2[z];
		}

		if (rdn5)// 
		{
			kid->traits.emig1.push_back(mummy.traits.emig1[z]); // 
			kid->traits.g_emig += mummy.traits.emig1[z]; // 

		}
		else {
			kid->traits.emig1.push_back(mummy.traits.emig2[z]);
			kid->traits.g_emig += mummy.traits.emig2[z];
		}


		// dads loci for emig
		if (rdn6)rdn6 -= recomb(rdgen); // 
		else rdn6 += recomb(rdgen);// 

		if (rdn6)//
		{
			kid->traits.emig2.push_back(daddy.traits.emig1[z]); // 
			kid->traits.g_emig += daddy.traits.emig1[z];

		}
		else {
			kid->traits.emig2.push_back(daddy.traits.emig2[z]);
			kid->traits.g_emig += daddy.traits.emig2[z];
		}

		// emigM

		if (rdn7)// mum loci
		{
			kid->traits.emigM1.push_back(mummy.traits.emigM1[z]); // 
			kid->traits.g_emigM += mummy.traits.emigM1[z]; // 

		}
		else {
			kid->traits.emigM1.push_back(mummy.traits.emigM2[z]);
			kid->traits.g_emigM += mummy.traits.emigM2[z];
		}


		// dads loci for emigM
		if (rdn8)rdn8 -= recomb(rdgen); // 
		else rdn8 += recomb(rdgen);// 

		if (rdn8)//
		{
			kid->traits.emigM2.push_back(daddy.traits.emigM1[z]); // 
			kid->traits.g_emigM += daddy.traits.emigM1[z];

		}
		else {
			kid->traits.emigM2.push_back(daddy.traits.emigM2[z]);
			kid->traits.g_emigM += daddy.traits.emigM2[z];
		}


	}

	// assign phenotype, which is equal to g 
	// pref

	kid->traits.p_pref = kid->traits.g_pref;

	// display
	if (kid->traits.g_display > 0)
	{
		kid->traits.p_display = kid->traits.g_display; // phenotype equals the genotype when trait is above 0
	}

	else
	{
		kid->traits.p_display = 0; // phenotype is bound by 0
	}

	// emig
	if (kid->traits.g_emig < 0)
	{
		kid->traits.p_emig = 0;
	}
	else 
	{
			if (kid->traits.g_emig > 0)
			{
			kid->traits.p_emig = 1.0;
			}
			else
			{
			kid->traits.p_emig = kid->traits.g_emig;
			}
	}

	// emigM

	if (kid->traits.g_emigM < 0)
	{
		kid->traits.p_emigM = 0;
	}
	else
	{
		if (kid->traits.g_emigM > 0)
		{
			kid->traits.p_emigM = 1.0;
		}
		else
		{
			kid->traits.p_emigM = kid->traits.g_emigM;
		}
	}


}

void mutation(int x,int y)// passing the information to access the population
{
	int Nmut;
	int ind;
	int locus;

	normal_distribution<> pref_mutdistr(0.0, (stdPref / sqrt(2.0*(double)Nloci))*0.1); //to scale the variation to the intial variation in the trait. The variation in the evvect size of the muation is 1/10 of the inital variation in the trait
	normal_distribution<> display_mutdistr(0.0, (stdDisplay / sqrt(2.0*(double)Nloci))*0.1);
	normal_distribution<> emig_mutdistr(0.0, (stdEmig / sqrt(2.0*(double)Nloci))*0.1);
	poisson_distribution<>mutdistr((double)popgrid[x][y].Noff*mutationrate*Nloci*2.0*4.0); //MT 02/02/21 : times 2 for diploid // additional because 3 traits
	uniform_int_distribution<>indistr(0, popgrid[x][y].Noff-1);
	uniform_int_distribution<>locidistr(0, Nloci * 8 - 1); // times six because 4 traits with each homolog
	
	Nmut = mutdistr(rdgen);


	//cout << "mutation starting" << endl;

	for (int z = 0; z < Nmut; z++)
	{
		ind = indistr(rdgen);// sampling the individual
		locus = locidistr(rdgen);// sampling which loci

		if (popgrid[x][y].Noff > 0)// added
		{


			if (locus < Nloci * 2) // deciding that its preference 
			{
				if (ind < popgrid[x][y].Foff)// deciding that its females with lower ind
				{

					popgrid[x][y].femaleOffspring[ind].mutation_effect(locus, Nloci, pref_mutdistr);
				}

				else
				{

					popgrid[x][y].maleOffspring[ind - popgrid[x][y].Foff].mutation_effect(locus, Nloci, pref_mutdistr);
				}
				//cout << "mutation effect Ok" << endl;
			}
			else
			{
				if (locus < Nloci * 4) // deciding that its display
				{
					if (ind < popgrid[x][y].Foff)// deciding that its females with lower ind
					{

						popgrid[x][y].femaleOffspring[ind].mutation_effect(locus, Nloci, display_mutdistr);
					}

					else
					{

						popgrid[x][y].maleOffspring[ind - popgrid[x][y].Foff].mutation_effect(locus, Nloci, display_mutdistr);
					}

				}
				else // emig
				{
					if (locus < Nloci * 6) // deciding that its emig
					{
						if (ind < popgrid[x][y].Foff)// deciding that its females with lower ind
						{

							popgrid[x][y].femaleOffspring[ind].mutation_effect(locus, Nloci, emig_mutdistr);
						}

						else
						{

							popgrid[x][y].maleOffspring[ind - popgrid[x][y].Foff].mutation_effect(locus, Nloci, emig_mutdistr);
						}

					}
					else // deciding emigM
					{
		
						if (ind < popgrid[x][y].Foff)// deciding that its females with lower ind
						{

							popgrid[x][y].femaleOffspring[ind].mutation_effect(locus, Nloci, emig_mutdistr);
						}

						else
						{

							popgrid[x][y].maleOffspring[ind - popgrid[x][y].Foff].mutation_effect(locus, Nloci, emig_mutdistr);
						}
					
					}
				}
				//cout << "mutation effect Ok" << endl;
			}
		}
	}

}

//--------------------------------SPACE-----------------------------------------------

void dispersal(void) {
	uniform_real_distribution<> Prob(0.0, 1.0);
	bernoulli_distribution costdistr(disp_cost);

	int new_x;
	int new_y;
	double distance, rnd_angle, R1, x_rand, y_rand;
	double emp; //emigration probabilty 

	vector<Individuals>::iterator iter;
	std::uniform_real_distribution<>unireal(0.000000001, 0.999);
	std::uniform_real_distribution<>position(0.0, 0.999); // to get the postiion of the ind in the cell

	
	//bernoulli_distribution dispdistr(disprob);

	

	for (int i = 0; i < xmax; ++i)
	{
		for (int j = 0; j < ymax; ++j) 
		{
			if (popgrid[i][j].N > 0)
			{
				//for female offspring
				for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++)
				{
				
					if(disp_evol) emp=iter->traits.p_emig;
					else emp= disprob;

					if(unireal(rdgen)<emp)// individual is dispersing
					{

						//make indivudal disperser type true

						iter->dispersal_type = true;

						if (!costdistr(rdgen)) // if ind suvives then. (Cost dispersal)
						{
							//cout << "female dispersing"<< iter->x << iter->y << endl;


							x_rand = position(rdgen);// sample random position
							y_rand = position(rdgen);
							do
							{

								do
								{
									R1 = unireal(rdgen);// random number
									distance = (-1.0*mean_distance)*log(R1);// calc distance
									rnd_angle = Prob(rdgen)*2.0*PI;// angle
									new_x = (int)(distance*cos(rnd_angle) / cell_resolution + x_rand + i); //translate the continouse distance into a new position in the discrete cell
									new_y = (int)(distance*sin(rnd_angle) / cell_resolution + y_rand + j);

								} while (new_x == i && new_y == j);

							} while (new_x<0 || new_x>xmax - 1 || new_y<0 || new_y>ymax - 1);

							// give indivudal newx the value of new_x

							iter->newx = new_x;
							iter->newy = new_y;

							if (dispersal_out)
							{


							

								if (g%out_ind_interval == 0)
								{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
								else
								{
									if (g > out2_limit && g%out_ind_interval2 == 0)
									{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
									}
								}

							}
							//cout << iter->x << iter->y << " and " << new_x << new_y << endl;
							iter->x = new_x;
							iter->y = new_y;

							popgrid[new_x][new_y].females.push_back(*iter);
							popgrid[new_x][new_y].N++;
							popgrid[i][j].N--;
							popgrid[new_x][new_y].Nf++;
							popgrid[i][j].Nf--;
						}
						else 
						{
							popgrid[i][j].N--; 
							popgrid[i][j].Nf--;
						}
					}
					else
					{

						if (dispersal_out)
						{

						
							//extract resident information
							iter->dispersal_type = false;
							iter->newx = iter->x;
							iter->newy = iter->y;

							if (g%out_ind_interval == 0)
							{
							iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
							else
							{
								if (g > out2_limit && g%out_ind_interval2 == 0)
								{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
							}

						}

						popgrid[i][j].females.push_back(*iter);
					}
				}
				// same for male offspring 
				for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++)
				{
					if (disp_evol)
					{
						if (disp_sex) emp = iter->traits.p_emigM; // male specific dispersal
						else emp = iter->traits.p_emig; // both sexes the same
					}
					else emp = disprob;

					if (unireal(rdgen) < emp)// individual is dispersing
					{
						// make ind a disperser type

						iter->dispersal_type = true;

						if (!costdistr(rdgen)) // if ind suvives then. (Cost dispersal)
						{
							//cout << "male disperser" << endl;


							x_rand = position(rdgen);// sample random position
							y_rand = position(rdgen);
							do
							{

								do
								{
									R1 = unireal(rdgen);// random number
									distance = (-1.0*mean_distance)*log(R1);// calc distance
									rnd_angle = Prob(rdgen)*2.0*PI;// angle
									new_x = (int)(distance*cos(rnd_angle) / cell_resolution + x_rand + i); //translate the continouse distance into a new position in the discrete cell
									new_y = (int)(distance*sin(rnd_angle) / cell_resolution + y_rand + j);

								} while (new_x == i && new_y == j);

							} while (new_x<0 || new_x>xmax - 1 || new_y<0 || new_y>ymax - 1);

							// give ind new coordinates 

							iter->newx = new_x;
							iter->newy = new_y;

							if (dispersal_out)
							{




								if (g%out_ind_interval == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
								else
								{
									if (g > out2_limit && g%out_ind_interval2 == 0)
									{
										iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
									}
								}

							}
							
							iter->x = new_x;
							iter->y = new_y;

							popgrid[new_x][new_y].males.push_back(*iter);
							popgrid[new_x][new_y].N++;
							popgrid[i][j].N--;
							popgrid[new_x][new_y].Nm++;
							popgrid[i][j].Nm--;
						}
						else
						{
							popgrid[i][j].N--;
							popgrid[i][j].Nm--;
						}
					}
					else
					{

						iter->dispersal_type = false;
						iter->newx = iter->x;
						iter->newy = iter->y;

						if (dispersal_out)
						{



							if (g%out_ind_interval == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
							else
							{
								if (g > out2_limit && g%out_ind_interval2 == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
							}
						}

						popgrid[i][j].males.push_back(*iter);
					}

				}


			
			}

			if(!popgrid[i][j].tmp_females.empty()) popgrid[i][j].tmp_females.clear();
			if (!popgrid[i][j].tmp_males.empty()) popgrid[i][j].tmp_males.clear();

			//cout << "nfemales " << popgrid[i][j].Nf << endl;
			//cout << "nmales " << popgrid[i][j].Nm << endl;


		}
	}
}
 // current dispersal function commented out

void dispersal_ds(void) {
	uniform_real_distribution<> Prob(0.0, 1.0);
	bernoulli_distribution costdistr(disp_cost);

	int new_x;
	int new_y;
	double distance, rnd_angle, R1, x_rand, y_rand;
	double emp; //emigration probabilty 

	vector<Individuals>::iterator iter;
	std::uniform_real_distribution<>unireal(0.000000001, 0.999);
	std::uniform_real_distribution<>position(0.0, 0.999); // to get the postiion of the ind in the cell


	//bernoulli_distribution dispdistr(disprob);



	for (int i = 0; i < xmax; ++i)
	{
		for (int j = 0; j < ymax; ++j)
		{
			if (popgrid[i][j].Noff > 0)
			{
				//for female offspring
				for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++)
				{

					if (disp_evol) emp = iter->traits.p_emig;
					else emp = disprob;

					if (unireal(rdgen) < emp)// individual is dispersing
					{

						//make indivudal disperser type true

						iter->dispersal_type = true;

						if (!costdistr(rdgen)) // if ind suvives then. (Cost dispersal)
						{
							//cout << "female dispersing"<< iter->x << iter->y << endl;


							x_rand = position(rdgen);// sample random position
							y_rand = position(rdgen);
							do
							{

								do
								{
									R1 = unireal(rdgen);// random number
									distance = (-1.0*mean_distance)*log(R1);// calc distance
									rnd_angle = Prob(rdgen)*2.0*PI;// angle
									new_x = (int)(distance*cos(rnd_angle) / cell_resolution + x_rand + i); //translate the continouse distance into a new position in the discrete cell
									new_y = (int)(distance*sin(rnd_angle) / cell_resolution + y_rand + j);

								} while (new_x == i && new_y == j);

							} while (new_x<0 || new_x>xmax - 1 || new_y<0 || new_y>ymax - 1);

							// give indivudal newx the value of new_x

							iter->newx = new_x;
							iter->newy = new_y;

							if (dispersal_out)
							{




								if (g%out_ind_interval == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
								else
								{
									if (g > out2_limit && g%out_ind_interval2 == 0)
									{
										iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
									}
								}

							}
							//cout << iter->x << iter->y << " and " << new_x << new_y << endl;
							iter->x = new_x;
							iter->y = new_y;

							popgrid[new_x][new_y].tmp_females.push_back(*iter);
							popgrid[new_x][new_y].Noff++;
							popgrid[i][j].Noff--;
							popgrid[new_x][new_y].Foff++;
							popgrid[i][j].Foff--;
						}
						else
						{
							popgrid[i][j].Noff--;
							popgrid[i][j].Foff--;
						}
					}
					else
					{

						if (dispersal_out)
						{


							//extract resident information
							iter->dispersal_type = false;
							iter->newx = iter->x;
							iter->newy = iter->y;

							if (g%out_ind_interval == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
							else
							{
								if (g > out2_limit && g%out_ind_interval2 == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
							}

						}

						popgrid[i][j].tmp_females.push_back(*iter);
					}
				}
				// same for male offspring 
				for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++)
				{
					
					if (disp_evol)
					{
						if (disp_sex) emp = iter->traits.p_emigM; // male specific dispersal
						else emp = iter->traits.p_emig; // both sexes the same
					}
					else emp = disprob;

					if (unireal(rdgen) < emp)// individual is dispersing
					{
						// make ind a disperser type

						iter->dispersal_type = true;

						if (!costdistr(rdgen)) // if ind suvives then. (Cost dispersal)
						{
							//cout << "male disperser" << endl;


							x_rand = position(rdgen);// sample random position
							y_rand = position(rdgen);
							do
							{

								do
								{
									R1 = unireal(rdgen);// random number
									distance = (-1.0*mean_distance)*log(R1);// calc distance
									rnd_angle = Prob(rdgen)*2.0*PI;// angle
									new_x = (int)(distance*cos(rnd_angle) / cell_resolution + x_rand + i); //translate the continouse distance into a new position in the discrete cell
									new_y = (int)(distance*sin(rnd_angle) / cell_resolution + y_rand + j);

								} while (new_x == i && new_y == j);

							} while (new_x<0 || new_x>xmax - 1 || new_y<0 || new_y>ymax - 1);

							// give ind new coordinates 

							iter->newx = new_x;
							iter->newy = new_y;

							if (dispersal_out)
							{




								if (g%out_ind_interval == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
								else
								{
									if (g > out2_limit && g%out_ind_interval2 == 0)
									{
										iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
									}
								}

							}

							iter->x = new_x;
							iter->y = new_y;

							popgrid[new_x][new_y].tmp_males.push_back(*iter);
							popgrid[new_x][new_y].Noff++;
							popgrid[i][j].Noff--;
							popgrid[new_x][new_y].Moff++;
							popgrid[i][j].Moff--;
						}
						else
						{
							popgrid[i][j].Noff--;
							popgrid[i][j].Moff--;
						}
					}
					else
					{

						iter->dispersal_type = false;
						iter->newx = iter->x;
						iter->newy = iter->y;

						if (dispersal_out)
						{



							if (g%out_ind_interval == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
							else
							{
								if (g > out2_limit && g%out_ind_interval2 == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
							}
						}

						popgrid[i][j].tmp_males.push_back(*iter);
					}

				}



			}

			if (!popgrid[i][j].femaleOffspring.empty()) popgrid[i][j].femaleOffspring.clear();
			if (!popgrid[i][j].maleOffspring.empty()) popgrid[i][j].maleOffspring.clear();

			//cout << "nfemales " << popgrid[i][j].Nf << endl;
			//cout << "nmales " << popgrid[i][j].Nm << endl;


		}
	}
}


void dispersal_control(void) {
	uniform_real_distribution<> Prob(0.0, 1.0);
	bernoulli_distribution costdistr(disp_cost);

	int new_x;
	int new_y;
	double distance, rnd_angle, R1, x_rand, y_rand;

	double sum_pref = 0; // sum of preference values for calculating allelic mean
	double sum_display = 0; // sum of display vlaues for calculating allelic mean

	double sd_pref = 0;
	double sd_display = 0;

	


	vector <int> sel_ind; // keeps ID of all individuals
	vector<Individuals>::iterator iter;
	std::uniform_real_distribution<>unireal(0.000000001, 1.0);
	std::uniform_real_distribution<>position(0.0, 0.999); // to get the postiion of the ind in the cell

	bernoulli_distribution dispdistr(disprob);

	for (int i = 0; i < xmax; ++i)
	{
		for (int j = 0; j < ymax; ++j)
		{

			//cout <<"NR F in pop vector" << popgrid[i][j].tmp_females.size() << endl;


			popgrid[i][j].bd_Nf = popgrid[i][j].Nf;
			popgrid[i][j].bd_Nm = popgrid[i][j].Nm;

		}

	}

	// looping thoru popgrid cells i and j, to calculate mean and sd of allelic deistribution
	for (int i = 0; i < xmax; ++i)
	{
		for (int j = 0; j < ymax; ++j)
		{
			if (popgrid[i][j].N > 0)
			{
				// calculate mean of populations allelic values by summing all the genotypic values and divide by number of inds*2*number of loci
				// then you have the distribution for xy that i can use


				// calc means
				for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++)
				{
					sum_pref += iter->traits.g_pref;
					sum_display += iter->traits.g_display;
				}

				for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++)
				{
					sum_pref += iter->traits.g_pref;
					sum_display += iter->traits.g_display;
				}

				//cout << " sum display : " << sum_display << " and sum pref: " << sum_pref << endl;
				
				popgrid[i][j].mean_pref = sum_pref / ((double)popgrid[i][j].N*2.0*(double)Nloci);
				popgrid[i][j].mean_display = sum_display / ((double)popgrid[i][j].N *2.0 *(double)Nloci);
				
				sum_pref = 0;
				sum_display = 0;
				
				//cout << "mean pref:" << popgrid[i][j].mean_pref << endl;
				//cout << "mean display:" << popgrid[i][j].mean_display << endl;


				//calc sds
				for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++)
				{
					sd_pref += pow(((iter->traits.g_pref/(2.0*(double)Nloci))-popgrid[i][j].mean_pref),2.0);
					sd_display += pow(((iter->traits.g_display/(2.0*(double)Nloci))- popgrid[i][j].mean_display),2.0);
					
				}

				for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++)
				{
					sd_pref += pow( ((iter->traits.g_pref/(2.0*(double)Nloci))- popgrid[i][j].mean_pref),2.0);
					sd_display += pow(((iter->traits.g_display/(2.0*(double)Nloci)) - popgrid[i][j].mean_display),2.0);
				}

				popgrid[i][j].sd_pref =sqrt(sd_pref/((double)popgrid[i][j].N-1.0));
				popgrid[i][j].sd_display = sqrt(sd_display/((double)popgrid[i][j].N-1.0));

				//cout << "sd pref:" << popgrid[i][j].sd_pref << endl;
				//cout << "sd display:" << popgrid[i][j].sd_display << endl;

				sd_display = 0;
				sd_pref = 0;
			}

		}
	}

	for (int i = 0; i < xmax; ++i)
	{
		for (int j = 0; j < ymax; ++j)
		{
			if (popgrid[i][j].N > 0)
			{



				//for female offspring
				for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++)
				{
					//cout << "tmpfemale pref value: " << iter->traits.p_pref << endl;

					if (dispdistr(rdgen))// individual is dispersing
					{

						//make indivudal disperser type true

						iter->dispersal_type = true;

						if (!costdistr(rdgen)) // if ind suvives then. (Cost dispersal)
						{
							//cout << "female dispersing"<< iter->x << iter->y << endl;


							x_rand = position(rdgen);// sample random position
							y_rand = position(rdgen);
							do
							{

								do
								{
									R1 = unireal(rdgen);// random number
									distance = (-1.0*mean_distance)*log(R1);// calc distance
									rnd_angle = Prob(rdgen)*2.0*PI;// angle
									new_x = (int)(distance*cos(rnd_angle) / cell_resolution + x_rand + i); //translate the continouse distance into a new position in the discrete cell
									new_y = (int)(distance*sin(rnd_angle) / cell_resolution + y_rand + j);

								} while (new_x == i && new_y == j);

							} while (new_x<0 || new_x>xmax - 1 || new_y<0 || new_y>ymax - 1);

							// give indivudal newx the value of new_x



							iter->newx = new_x;
							iter->newy = new_y;


							if (g%out_ind_interval == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
							else
							{
								if (g > out2_limit && g%out_ind_interval2 == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
							}


							//cout << " Old coord xy "<< iter->x << iter->y << " and new  " << new_x << new_y << endl;
							iter->x = new_x;
							iter->y = new_y;

							if (popgrid[new_x][new_y].N > 0) // change genetic and phenotypic values of individual by sampling from new residents, when there are residents available in the before dispersal populations
							{
								


								if (!iter->traits.pref1.empty()) iter->traits.pref1.clear();// clear the chromosomes of the individual
								if (!iter->traits.pref2.empty()) iter->traits.pref2.clear();
								if (!iter->traits.display1.empty()) iter->traits.display1.clear();
								if (!iter->traits.display2.empty()) iter->traits.display2.clear();

								if (popgrid[new_x][new_y].sd_display > 0) // for cases where sd in display is above zero
								{
									if (popgrid[new_x][new_y].sd_pref > 0)
									{

										normal_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].sd_display);
										normal_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].sd_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}

										

										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}

									else
									{ 
										normal_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].sd_display);
										uniform_real_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].mean_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}


										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}



								}

								else // if sd_display 0
								{

									if (popgrid[new_x][new_y].sd_pref > 0)
									{

										uniform_real_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].mean_display);
										normal_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].sd_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}


										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}


									else
									{
										uniform_real_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].mean_display);
										uniform_real_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].mean_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}


										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}
								}
							}


							else
							{
								
								//cout << indID <<"NO IND THERE" << endl;
							}

							popgrid[new_x][new_y].females.push_back(*iter);
							popgrid[new_x][new_y].N++;
							popgrid[i][j].N--;
							popgrid[new_x][new_y].Nf++;
							popgrid[i][j].Nf--;

						}
						else
						{
							popgrid[i][j].N--;
							popgrid[i][j].Nf--;
						}
					}

					else
					{
						//extract resident information
						iter->dispersal_type = false;
						iter->newx = iter->x;
						iter->newy = iter->y;

						if (g%out_ind_interval == 0)
						{
							iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
						}
						else
						{
							if (g > out2_limit && g%out_ind_interval2 == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
						}


						popgrid[i][j].females.push_back(*iter);



					}
				}
			
				// same for male offspring 
				for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++)
				{
					if (dispdistr(rdgen))// individual is dispersing
					{
						// make ind a disperser type

						iter->dispersal_type = true;

						if (!costdistr(rdgen)) // if ind suvives then. (Cost dispersal)
						{
							//cout << "male disperser" << endl;


							x_rand = position(rdgen);// sample random position
							y_rand = position(rdgen);
							do
							{

								do
								{
									R1 = unireal(rdgen);// random number
									distance = (-1.0*mean_distance)*log(R1);// calc distance
									rnd_angle = Prob(rdgen)*2.0*PI;// angle
									new_x = (int)(distance*cos(rnd_angle) / cell_resolution + x_rand + i); //translate the continouse distance into a new position in the discrete cell
									new_y = (int)(distance*sin(rnd_angle) / cell_resolution + y_rand + j);

								} while (new_x == i && new_y == j);

							} while (new_x<0 || new_x>xmax - 1 || new_y<0 || new_y>ymax - 1);

							// give ind new coordinates 

							iter->newx = new_x;
							iter->newy = new_y;


							if (g%out_ind_interval == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
							else
							{
								if (g > out2_limit && g%out_ind_interval2 == 0)
								{
									iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
								}
							}


							iter->x = new_x;
							iter->y = new_y;

							// transfer migrants genetic and phenotypic values with randomly selected resident 
							//cout << popgrid[new_x][new_y].bd_Nm << " Nms " << endl; //checked

							// if the new population where indivudal is going has some inidivuals in there, sample from its allelic distribution

							if (popgrid[new_x][new_y].N > 0) // change genetic and phenotypic values of individual by sampling from new residents, when there are residents available in the before dispersal populations
							{

								

								if (!iter->traits.pref1.empty()) iter->traits.pref1.clear();// clear the chromosomes of the individual
								if (!iter->traits.pref2.empty()) iter->traits.pref2.clear();
								if (!iter->traits.display1.empty()) iter->traits.display1.clear();
								if (!iter->traits.display2.empty()) iter->traits.display2.clear();

								if (popgrid[new_x][new_y].sd_display > 0) // for cases where sd in display is above zero
								{
									if (popgrid[new_x][new_y].sd_pref > 0)
									{

										normal_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].sd_display);
										normal_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].sd_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}

										

										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}

									else
									{// in case female allelic sd is 0
										normal_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].sd_display);
										uniform_real_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].mean_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}


										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}



								}

								else // if sd_display 0
								{

									if (popgrid[new_x][new_y].sd_pref > 0)
									{

										uniform_real_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].mean_display);
										normal_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].sd_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}


										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}


									else
									{// in case female allelic sd is 0 adn male allelic display is 0
										uniform_real_distribution <> display_alleledistr(popgrid[new_x][new_y].mean_display, popgrid[new_x][new_y].mean_display);
										uniform_real_distribution <> pref_alleledistr(popgrid[new_x][new_y].mean_pref, popgrid[new_x][new_y].mean_pref);

										for (int i = 0; i < Nloci; i++)
										{
											iter->traits.display1.push_back(display_alleledistr(rdgen)); // fill chromosome up with alleles for loci drwan from a random gen
											iter->traits.display2.push_back(display_alleledistr(rdgen));

											iter->traits.pref1.push_back(pref_alleledistr(rdgen));
											iter->traits.pref2.push_back(pref_alleledistr(rdgen));

											iter->traits.g_display += iter->traits.display1[i] + iter->traits.display2[i]; // create additive trait value with equal effects of both chromosomes
											iter->traits.g_pref += iter->traits.pref1[i] + iter->traits.pref2[i];

										}


										if (iter->traits.g_display > 0)
										{
											iter->traits.p_display = iter->traits.g_display; // phenotype equals the genotype when trait is above 0
										}

										else
										{
											iter->traits.p_display = 0; // phenotype is bound by 0
										}

										iter->traits.p_pref = iter->traits.g_pref;

									}

								}
							}
						
					

							else
							{
							
							//cout << indID << "NO IND THERE" << endl;
							}


							popgrid[new_x][new_y].males.push_back(*iter);
							popgrid[new_x][new_y].N++;
							popgrid[i][j].N--;
							popgrid[new_x][new_y].Nm++;
							popgrid[i][j].Nm--;
						}
						else
						{
							popgrid[i][j].N--;
							popgrid[i][j].Nm--;
						}
					}
					else
					{

						iter->dispersal_type = false;
						iter->newx = iter->x;
						iter->newy = iter->y;

						if (g%out_ind_interval == 0)
						{
							iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
						}
						else
						{
							if (g > out2_limit && g%out_ind_interval2 == 0)
							{
								iter->outind_dispersal(r, g, &inds_dispersal); // output  large format
							}
						}
						popgrid[i][j].males.push_back(*iter);
					}

				}

			}

			//if (!popgrid[i][j].tmp_females.empty()) popgrid[i][j].tmp_females.clear();
			//if (!popgrid[i][j].tmp_males.empty()) popgrid[i][j].tmp_males.clear();

			//cout << "nfemales " << popgrid[i][j].Nf << endl;
			//cout << "nmales " << popgrid[i][j].Nm << endl;


		}
	}

	for (int i = 0; i < xmax; ++i)
	{
		for (int j = 0; j < ymax; ++j)
		{
			if (popgrid[i][j].N > 0)
			{

				if (!popgrid[i][j].tmp_females.empty()) popgrid[i][j].tmp_females.clear();
				if (!popgrid[i][j].tmp_males.empty()) popgrid[i][j].tmp_males.clear();

			}
		}
	}
}


void survival(void) {
	vector<Individuals>::iterator iter;
	
	double sp; // survival probability

	double sum_vf = 0.0; // sum of all the females viablity in a population
	double sum_vm = 0.0; // sum of all the males viablity in a population

	for (int i = 0; i < xmax; i++)
	{
		for (int j = 0; j < ymax; j++)
		{
			sum_vf = 0.0;
			sum_vm = 0.0;

			if (popgrid[i][j].Noff > 0)
			{

				//TRAIT DEPENDENT COST (get indi into temp_vector 2)




				// FEMALE (first)
				
				if (f_costs) // if female pref is costly, and hence selection is acting on the trait
				{
					if (fsel_din) // selection is density-independent
					{



						if (fc_hard)// and if selection is hard
						{
							for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0)));
								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;
								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}

							}

						}

						else //and  if selection is soft ( and density-independent), has to be worked on!!!!!!!!!!!!!!!!!!
						{
							for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0)));
								sum_vf += iter->traits.viability;

								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;
							}

							for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++)
							{
								iter->traits.viability = (iter->traits.viability / sum_vf);
								bernoulli_distribution selection_dist(iter->traits.viability);
								//cout << "female trait2 " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;

								if (selection_dist(rdgen))
								{

									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}
							}
						}

					}

					else // if selection is density-dependent then
					{
						if (fc_hard)// and if selection is hard
						{
							for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability =fmin((grid[i][j].K/(double)popgrid[i][j].Noff)* exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0))),1.0); // having viability multiplied by density dependence
								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;
								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}

							}

						}
						else // and if selection is soft (density-dependent)
						{
							//calculate the viability of each female
							for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = fmin(exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0))), 1.0); // having viability multiplied by density dependence
								sum_vf += iter->traits.viability;
								
								//cout << sum_vf << endl;
								//cout << iter->traits.viability << endl;
								//cout << "end" << endl;

								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;

							}

							for(iter=popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end();iter++)
							{
								
								sp = fmin((iter->traits.viability / sum_vf)*grid[i][j].K,1);

								//cout <<"prob=" << sp << endl;
								bernoulli_distribution selection_dist(sp);
							

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}

							}

						}





					}
				}

				else // if no costs apply to females
				{

					for (iter = popgrid[i][j].femaleOffspring.begin(); iter != popgrid[i][j].femaleOffspring.end(); iter++) // calculate viability for each female offspring 
					{
						iter->traits.viability = 1.0;

						popgrid[i][j].tmp2_females.push_back(*iter);
						popgrid[i][j].Noff2++;
						popgrid[i][j].Foff2++;

					}
				}
					
				// MALES  (second)

				if (m_costs)
				{

					if (msel_din)// if selection on male is density-independent
					{



						if (mc_hard)
						{
							for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++) // calculate viability for each male offspring 
							{
								iter->traits.viability =fmin(exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0))),1);

								//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;

								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}
							}
						}

						else // selection is soft. NEEDS WORK !!!!!!!!!!!!!!!!!!!!!!!
						{
							for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++) // calculate viability for each male offspring 
							{
								iter->traits.viability = exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0)));
								sum_vm += iter->traits.viability;

								//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;
							}


							//cout << " males selection weak: sum_vm " << sum_vm << endl;

							for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++)
							{
								iter->traits.viability = fmin((iter->traits.viability / sum_vm), 1);
								bernoulli_distribution selection_dist(iter->traits.viability);

								//cout << "viablity: " << iter->traits.viability << endl;

								if (selection_dist(rdgen))
								{

									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}

							}
						}
					}
					
					else // if selection in density-dependent
					{ 
						if (mc_hard) // desnity dependent hard selection
						{

							for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++) // calculate viability for each male offspring 
							{
								iter->traits.viability =fmin((grid[i][j].K /(double)popgrid[i][j].Noff)*(exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0)))),1.0);

								//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;

								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}
							}
						}

						else // if selection is density dependent soft
						{
							//calculate the viability of each male
							for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = fmin(exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0))), 1.0); // having viability multiplied by density dependence
								sum_vm += iter->traits.viability;
								//cout << sum_vm << endl;
								//cout << iter->traits.viability << endl;
								//cout << "end" << endl;


							}

							for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++)
							{

								sp = fmin((iter->traits.viability / sum_vm)*grid[i][j].K, 1);

								//cout << "prob=" << sp << endl;
								bernoulli_distribution selection_dist(sp);


								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}

							}

						}
					}
				}

				else // if no costs apply to males
				{

					for (iter = popgrid[i][j].maleOffspring.begin(); iter != popgrid[i][j].maleOffspring.end(); iter++) // calculate viability for each female offspring 
					{
						iter->traits.viability = 1.0;

						//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;
						
						popgrid[i][j].tmp2_males.push_back(*iter);
						popgrid[i][j].Noff2++;
						popgrid[i][j].Moff2++;

					}
				}

				//cout << popgrid[i][j].Noff2 << endl;
				//cout << popgrid[i][j].Moff2 << endl;
				//cout << popgrid[i][j].Foff2 << endl;

				// TRAIT DEPENEDNT SELECTION OVER. Moving ind from tmp_2 to tmp_ vector

				sp = fmin(grid[i][j].K / (double)popgrid[i][j].Noff2, 1.0); // survival probability based on the new offspring vector, hence, prior mortalities make room for the rest
				bernoulli_distribution surv_distr(sp);

				//FEMALES FIRST

				
					
					for (iter = popgrid[i][j].tmp2_females.begin(); iter != popgrid[i][j].tmp2_females.end(); iter++) // for females in the population
					{
						//cout << "SURB+VIVALfemale pref value: " << iter->traits.p_pref << endl;
						if (surv_distr(rdgen))
						{
							popgrid[i][j].tmp_females.push_back(*iter);
							popgrid[i][j].Nf++;
							popgrid[i][j].N++;
						}
					}

				
				

			

				
				
				// MALES

			
				
					for (iter = popgrid[i][j].tmp2_males.begin(); iter != popgrid[i][j].tmp2_males.end(); iter++) // for males in the population
					{
						//cout << "male pref value: " << iter->traits.p_display << endl;
						if (surv_distr(rdgen))
						{
							popgrid[i][j].tmp_males.push_back(*iter);
							popgrid[i][j].Nm++;
							popgrid[i][j].N++;
						}
					}

				



			}

			if (!popgrid[i][j].femaleOffspring.empty()) popgrid[i][j].femaleOffspring.clear();// clear temporary vectors
			if (!popgrid[i][j].tmp2_females.empty()) popgrid[i][j].tmp2_females.clear();
			if (!popgrid[i][j].maleOffspring.empty()) popgrid[i][j].maleOffspring.clear();
			if (!popgrid[i][j].tmp2_males.empty()) popgrid[i][j].tmp2_males.clear();
			


			popgrid[i][j].Noff = 0;
			popgrid[i][j].Foff = 0;
			popgrid[i][j].Moff = 0;
			popgrid[i][j].Noff2 = 0;
			popgrid[i][j].Foff2 = 0;
			popgrid[i][j].Moff2 = 0;
			
		}
	}
}

void survival_ds(void) {
	vector<Individuals>::iterator iter;

	double sp; // survival probability

	double sum_vf = 0.0; // sum of all the females viablity in a population
	double sum_vm = 0.0; // sum of all the males viablity in a population

	for (int i = 0; i < xmax; i++)
	{
		for (int j = 0; j < ymax; j++)
		{
			sum_vf = 0.0;
			sum_vm = 0.0;

			if (popgrid[i][j].Noff > 0)
			{

				//TRAIT DEPENDENT COST (get indi into temp_vector 2)




				// FEMALE (first)

				if (f_costs) // if female pref is costly, and hence selection is acting on the trait
				{
					if (fsel_din) // selection is density-independent
					{



						if (fc_hard)// and if selection is hard
						{
							for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0)));
								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;
								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives F" << endl;
									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}

							}

						}

						else //and  if selection is soft ( and density-independent), has to be worked on!!!!!!!!!!!!!!!!!!
						{
							for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0)));
								sum_vf += iter->traits.viability;

								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;
							}

							for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++)
							{
								iter->traits.viability = (iter->traits.viability / sum_vf);
								bernoulli_distribution selection_dist(iter->traits.viability);
								//cout << "female trait2 " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;

								if (selection_dist(rdgen))
								{

									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}
							}
						}

					}

					else // if selection is density-dependent then
					{
						if (fc_hard)// and if selection is hard
						{
							for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = fmin((grid[i][j].K / (double)popgrid[i][j].Noff)* exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0))), 1.0); // having viability multiplied by density dependence
								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;
								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}

							}

						}
						else // and if selection is soft (density-dependent)
						{
							//calculate the viability of each female
							for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = fmin(exp((-1.0*std::pow(iter->traits.g_pref - optima_f, 2.0)) / (2 * std::pow(w_f, 2.0))), 1.0); // having viability multiplied by density dependence
								sum_vf += iter->traits.viability;

								//cout << sum_vf << endl;
								//cout << iter->traits.viability << endl;
								//cout << "end" << endl;

								//cout << "female trait " << iter->traits.g_pref << "   viablity: " << iter->traits.viability << endl;

							}

							for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++)
							{

								sp = fmin((iter->traits.viability / sum_vf)*grid[i][j].K, 1);

								//cout <<"prob=" << sp << endl;
								bernoulli_distribution selection_dist(sp);


								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_females.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Foff2++;

								}

							}

						}





					}
				}

				else // if no costs apply to females
				{

					for (iter = popgrid[i][j].tmp_females.begin(); iter != popgrid[i][j].tmp_females.end(); iter++) // calculate viability for each female offspring 
					{
						iter->traits.viability = 1.0;

						popgrid[i][j].tmp2_females.push_back(*iter);
						popgrid[i][j].Noff2++;
						popgrid[i][j].Foff2++;

					}
				}

				// MALES  (second)

				if (m_costs)
				{

					if (msel_din)// if selection on male is density-independent
					{



						if (mc_hard)
						{
							for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++) // calculate viability for each male offspring 
							{
								iter->traits.viability = fmin(exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0))), 1);

								//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;

								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives M" << endl;
									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}
							}
						}

						else // selection is soft. NEEDS WORK !!!!!!!!!!!!!!!!!!!!!!!
						{
							for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++) // calculate viability for each male offspring 
							{
								iter->traits.viability = exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0)));
								sum_vm += iter->traits.viability;

								//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;
							}


							//cout << " males selection weak: sum_vm " << sum_vm << endl;

							for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++)
							{
								iter->traits.viability = fmin((iter->traits.viability / sum_vm), 1);
								bernoulli_distribution selection_dist(iter->traits.viability);

								//cout << "viablity: " << iter->traits.viability << endl;

								if (selection_dist(rdgen))
								{

									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}

							}
						}
					}

					else // if selection in density-dependent
					{
						if (mc_hard) // desnity dependent hard selection
						{

							for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++) // calculate viability for each male offspring 
							{
								iter->traits.viability = fmin((grid[i][j].K / (double)popgrid[i][j].Noff)*(exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0)))), 1.0);

								//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;

								bernoulli_distribution selection_dist(iter->traits.viability);

								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}
							}
						}

						else // if selection is density dependent soft
						{
							//calculate the viability of each male
							for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++) // calculate viability for each female offspring 
							{
								iter->traits.viability = fmin(exp((-1.0*std::pow(iter->traits.g_display - optima_m, 2.0)) / (2 * std::pow(w_m, 2.0))), 1.0); // having viability multiplied by density dependence
								sum_vm += iter->traits.viability;
								//cout << sum_vm << endl;
								//cout << iter->traits.viability << endl;
								//cout << "end" << endl;


							}

							for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++)
							{

								sp = fmin((iter->traits.viability / sum_vm)*grid[i][j].K, 1);

								//cout << "prob=" << sp << endl;
								bernoulli_distribution selection_dist(sp);


								if (selection_dist(rdgen))
								{
									//cout << "survives" << endl;
									popgrid[i][j].tmp2_males.push_back(*iter);
									popgrid[i][j].Noff2++;
									popgrid[i][j].Moff2++;

								}

							}

						}
					}
				}

				else // if no costs apply to males
				{

					for (iter = popgrid[i][j].tmp_males.begin(); iter != popgrid[i][j].tmp_males.end(); iter++) // calculate viability for each female offspring 
					{
						iter->traits.viability = 1.0;

						//cout << "male trait " << iter->traits.g_display << "   viablity: " << iter->traits.viability << endl;

						popgrid[i][j].tmp2_males.push_back(*iter);
						popgrid[i][j].Noff2++;
						popgrid[i][j].Moff2++;

					}
				}

				//cout <<"N2 "<< popgrid[i][j].Noff2 << endl;
				//cout << "moff2" << popgrid[i][j].Moff2 << endl;
				//cout << "foff2" << popgrid[i][j].Foff2 << endl;

				// TRAIT DEPENEDNT SELECTION OVER. Moving ind from tmp_2 to adult vector

				sp = fmin(grid[i][j].K / (double)popgrid[i][j].Noff2, 1.0); // survival probability based on the new offspring vector, hence, prior mortalities make room for the rest
				bernoulli_distribution surv_distr(sp);
				//cout << "pop reg: sp =" << sp << endl;

				//FEMALES FIRST



				for (iter = popgrid[i][j].tmp2_females.begin(); iter != popgrid[i][j].tmp2_females.end(); iter++) // for females in the population
				{
					//cout << "SURB+VIVALfemale pref value: " << iter->traits.p_pref << endl;
					if (surv_distr(rdgen))
					{
						popgrid[i][j].females.push_back(*iter);
						popgrid[i][j].Nf++; // 
						popgrid[i][j].N++;
					}
				}

				// MALES

				for (iter = popgrid[i][j].tmp2_males.begin(); iter != popgrid[i][j].tmp2_males.end(); iter++) // for males in the population
				{
					//cout << "male pref value: " << iter->traits.p_display << endl;
					if (surv_distr(rdgen))
					{
						popgrid[i][j].males.push_back(*iter);
						popgrid[i][j].Nm++;
						popgrid[i][j].N++;
					}
				}

			}

			if (!popgrid[i][j].tmp_females.empty()) popgrid[i][j].tmp_females.clear();// clear temporary vectors
			if (!popgrid[i][j].tmp2_females.empty()) popgrid[i][j].tmp2_females.clear();
			if (!popgrid[i][j].tmp_males.empty()) popgrid[i][j].tmp_males.clear();
			if (!popgrid[i][j].tmp2_males.empty()) popgrid[i][j].tmp2_males.clear();



			popgrid[i][j].Noff = 0;
			popgrid[i][j].Foff = 0;
			popgrid[i][j].Moff = 0;
			popgrid[i][j].Noff2 = 0;
			popgrid[i][j].Foff2 = 0;
			popgrid[i][j].Moff2 = 0;

		}
	}
}
//-------------------------------Output------------------------------------------------
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x))return"ERROR";
	return o.str();
}

void out_pop_header(void) {
	string name;

	name = dirout + "Sim" + Int2Str(simNr) + "_Pops.txt";
	pops.open(name.c_str());

	pops << "rep\tgen\tcoordx\tcoordy\tK\tS\tNf\tNm" << endl;

	
}

void out_ind_header(void)
{
	string name;
	
	name = dirout + "Sim" + Int2Str(simNr) + "_Inds.txt";
	inds.open(name.c_str());

	if (indout_slim)
	{
		inds << "rep\tgen\tcoordx\tcoordy\tg_pref\tg_display\tg_emig\tg_emigM\tsex" << endl;
	}

	else
	{
		inds << "rep\tgen\tcoordx\tcoordy\tsex\tg_pref\tp_pref\tg_display\tp_display\tg_emig\tp_emig\tg_emigM\tp_emigM" << endl;
	}
}

void out_inddispersal_header(void)
{
	string name;

	name = dirout + "Sim" + Int2Str(simNr) + "_Inds_dispersal.txt";
	inds_dispersal.open(name.c_str());

	inds_dispersal << "rep\tgen\tcoordx\tcoordy\tnewx\tnewy\tsex\tdispersal_status\tg_pref\tg_display" << endl;
}


void out_param_header(void) {
	string name;

	name = dirout + "Sim" + Int2Str(simNr) + "_parameters.txt";
	param.open(name.c_str());

	param << "SimNr\trep\tgen\txmax\tymax\tKmin\tKmax\tSmin\tSmax\thabturn\tdisprob\tmdistance\tdispcost\tdispevol\tdispsex\tNloci\tmeanPref\tstdPref\tmeanDisplay\tstdDisplay\tmeanEmig\tstdEmig\tmeanEmigM\tstdEmigM\tmutrate\tsampleallmates\tind_sample\tquantileS\tsubsamplemates\toffspring\toutpop\toutind\tcost_f\tfc_hard\tfsel_din\tw_f\topt_f\tcost_m\tmc_hard\tmsel_din\tw_m\topt_m" << endl;
}

void outparam(int simNr,std::ofstream*out) {

	*out << simNr << "\t" << rep << "\t" << gen << "\t" << xmax << "\t" << ymax << "\t"<< Kmin << "\t" << Kmax<< "\t"<< Smin <<"\t"<< Smax << "\t" << hab_turn << "\t" << disprob << "\t" << mean_distance << "\t" << disp_cost <<  "\t" << disp_evol << "\t" << disp_sex << "\t" << Nloci << "\t" << meanPref << "\t" << stdPref << "\t" << meanDisplay << "\t" << stdDisplay << "\t" <<meanEmig << "\t"<< stdEmig<< "\t" << meanEmigM << "\t" << stdEmigM << "\t" << mutationrate << "\t" << sample_all_mates << "\t" << ind_sample<<"\t"<< quantileS << "\t" << subsample_mates << "\t" << offspring << "\t" << out_pop_interval << "\t" << out_ind_interval <<"\t"<< f_costs<< "\t" << fc_hard << "\t" << fsel_din << "\t"<<w_f<<"\t"<<optima_f<<"\t"<<m_costs<<"\t"<< mc_hard <<"\t"<< msel_din<<"\t"<<w_m<<"\t"<<optima_m<< endl;
}

void sim_out(void) {
	
	vector<Individuals>::iterator iter;

		for (int i = 0; i < xmax; i++)
		{
			for (int j = 0; j < ymax; j++)
			{
					
					for (iter = popgrid[i][j].females.begin(); iter != popgrid[i][j].females.end(); iter++) // for females in the population
					{
						if (g < out2_limit && g%out_ind_interval == 0)
						{
							if (indout_slim) //output females
							{
								iter->outind_slim(r, g, &inds); // output small format
							}
							
							else
							{
								iter->outind(r, g, &inds); // output  large format
							}
							 
						}
						
						else
						{
							if (g > out2_limit && g%out_ind_interval2 == 0)
							{
								if (indout_slim) 
								{
									iter->outind_slim(r, g, &inds); // output small format
								}

								else
								{
									iter->outind(r, g, &inds); // output  large format
								}
							}
						}
					}

					for (iter = popgrid[i][j].males.begin(); iter != popgrid[i][j].males.end(); iter++) // for males in the population
					{
						if (g< out2_limit && g%out_ind_interval == 0)
						{
							if (indout_slim) //output males
							{
								iter->outind_slim(r, g, &inds); // output small format
							}

							else
							{
								iter->outind(r, g, &inds); // output  large format
							}
						}

						else
						{
							if (g > out2_limit && g%out_ind_interval2 == 0)
							{
								if (indout_slim) //
								{
									iter->outind_slim(r, g, &inds); // output small format
								}

								else
								{
									iter->outind(r, g, &inds); // output  large format
								}
							}
						}



					}
			
				
				if (g%out_pop_interval == 0 || (g > out2_limit && g%out_pop_interval2 == 0)) popgrid[i][j].outpop(r, g, i, j, grid[i][j].K,grid[i][j].S, &pops); // output population

			}
		}
	
}
