#include <iostream>
#include "../libs/QuantumParticle.h"

int main(int argc, char **argv)
{
	QuantumParticle::Particle<double> p;
	p.InitWaveFunction(1);
	// p.FlushToFiles();
	p.Simulation(0.0001,10,{1,4,9});	
	p.StabilityTest(1, std::pow(10,-7), std::pow(10,-3), 100, 0.01);
	return 1;
}
