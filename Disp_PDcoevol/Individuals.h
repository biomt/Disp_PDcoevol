#pragma once
// include library before internal scripts

#include<random>
#include<iostream>
#include<stdio.h>
#include<fstream>

//  headers

#include "traits.h"
using namespace std;



class Individuals {
public:
	Individuals(bool,int,int);
	~Individuals();

	bool sex;
	bool dispersal_type;
	int newx, newy;
	int x, y; // this are the coordinates for the ind 
	trait traits;// create an instance of the trait value 

	void setgenes(int L, std::normal_distribution<>NormPref, std::normal_distribution<>NormDisplay, std::normal_distribution<>NormEmig, std::normal_distribution<>NormEmigM);
	void mutation_effect(int L, int Nloci, normal_distribution<>effect);
	void outind(int repl, int gen, std::ofstream *out);
	void outind_slim(int repl, int gen, std::ofstream *out); // less output option to keep size of ind files down (no sex, no phenotype)
	void outind_dispersal(int repl, int gen, std::ofstream *out);
};
