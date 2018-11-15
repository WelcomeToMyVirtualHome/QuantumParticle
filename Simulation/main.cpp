#include <iostream>
#include "../libs/QuantumParticle.h"

int main(int argc, char **argv)
{
	QuantumParticle::Particle<double> p;
	p.InitWaveFunction(4);
	p.FlushToFiles();
	p.Simulation(0.0001,1,{1,4,9});
	// p.StabilityTest(1, std::pow(10,-8), std::pow(10,-4), 100, 0.01);
	return 1;
}
