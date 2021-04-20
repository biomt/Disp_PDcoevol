#pragma once

#define CLUSTER 0

#include <stdio.h>
#include <stdlib.h>

#if CLUSTER 

#include <unistd.h>

#else

#include <tchar.h>
#include <direct.h>
#include <io.h>

#endif

#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>

//include your own Header files and classes

#include"Landscape.h"
#include "Population.h"


using namespace std;

// need to have PI there for spatial part

const double PI = 3.141592654;

// global variables (parameters) 

//
int simNr =10; // change that to whaever simulation you are running, otherwise it will overwrite output!!!!!!!!!!!!

int rep = 2; // replicates of your simulations
int gen = 10; // number of generations

int r, g; // current replicate and current generation in the simulation

//Landscape

  const int xmax = 2;
  const int ymax = 2;

 const double cell_resolution = 100.0; // resolution in the cell

double hab_turn = 0.0;// rate at which habitat is destroyed (*indivuduals are killed off in sub population).



// Sampling of males by females


bool sample_all_mates = false; // if true sample all males, if false, subsample 
bool ind_sample = true; // if subsampling should be done in number of males (true) or otherwise use proportions to caclulate number
bool quantileS = false;	// if percentage is subsampled, (if true) then percentagage sampled in steps of 10, also S distrubution has to be from 1 to 10
int subsample_mates =1; // number of males being subsampled


// Reproduction, Mutatuion

double offspring=4.0; // mean offspring

double recprob = 0.5; // probabilty of recombination
double mutationrate = 0.001; // 

// Dispersal 


double disprob = 0.0001;// probability for dispersal deterministic
double mean_distance = 100.0; // mean dispersal distance in continouse space
double disp_cost = 0.0; // cost of dispersal

bool disp_evol = false; // if true, dispersal probability can evolve
bool disp_sex = false; // if true, sex-specific dispersal (males express EmigM)
bool dispersal_cntrl= false; // if true, dispersal_control function is run, otherwise the normal dispersal function is run
bool d_post = false; // post- survival dispersal

// Trait and Selection
int Nloci = 10;
double meanPref = 0.0;
double meanDisplay = 0.0;
double meanEmig = 0.5;
double meanEmigM = 0.5;


double stdPref = 1.0;
double stdDisplay = 1.0;
double stdEmig = 1.0;
double stdEmigM = 1.0;


bool m_costs = true; // if male trait is costly or not
bool mc_hard = true; // if selection on male trait is hard or soft
double w_m = 10.0; // strength of natural selection on male trait- will be squared 
double optima_m = 0.0;// global optima for male trait (could in the future be changed into a property of grid)

bool f_costs = true; // if female pref is costly or not
bool fc_hard = true; // if selection on female trait is hard or soft
double w_f = 10.0; // strength o natural selection on female trait
double optima_f = 0.0; // global optima for female pref


bool msel_din = true; // false if density-dependent selection on males occured
bool fsel_din = true; //false if density-dependent selection on females occured


string dir, dirout;
clock_t extime;
Landscape grid[xmax][ymax]; // Array with two dimensions for the space. using xmax or ymax to have these values clearly defined as constants.
Population popgrid[xmax][ymax]; // contains the populations

// declare functions

void runmodel(void);
void start(void);
void reproduction2(void);
void inheritance(Individuals*,Individuals,Individuals);// offspring and parents of calss individual go into the inheritance function
void mutation(int,int);
void dispersal(void);
void dispersal_ds(void);
void dispersal_control(void);
void survival(void);
void survival_ds(void);
void out_pop_header(void);// creates file, opens the file and writes the names of the first column (the headers)
void out_ind_header(void);
void out_inddispersal_header(void); // file
const string Int2Str(const int x);// to make a int into string
void outparam(int simNr,std::ofstream*out);
void out_param_header(void);
void sim_out(void);

// Random number generator

random_device rd;
mt19937 rdgen(rd()); // random number gnenerator using the rd. look for library random

//outputs

ofstream inds, inds_dispersal, pops, param;
int out_pop_interval = 20;
int out_ind_interval = 20;

int out2_limit = 80;

int out_pop_interval2 = 20;
int out_ind_interval2 = 20;


bool indout_slim = true;
bool dispersal_out = false;

// parameters for distributions

int Kmin = 200;
int Kmax = 200;

int Smin = 1; // if quantile is true : values from 1-10, if false: values from 1-100
int Smax = 1;



