#pragma once

class Landscape {

public:
	Landscape(); // constructor
	~Landscape(); // Deconstructor
	bool suitable; // contains suitable habitat that can contain a population
	double K; // carrying cap
	double S; // percentage of subsampled males in the population
};

