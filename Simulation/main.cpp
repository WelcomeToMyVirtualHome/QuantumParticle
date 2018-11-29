#include <iostream>
#include "../libs/QuantumParticle.h"

int main(int argc, char **argv)
{
	QuantumParticle::Particle<double> p;
	// p.SimulationOmegas(0.00008, 20, {1,4,9});
	// p.StabilityTest(1, std::pow(10,-8), std::pow(10,-4) * 1.2, 20);
	// p.Simulation(0.00002,20,{1},3*M_PI*M_PI/2);
	p.SimulationResonance(0.00002,20,1,3*M_PI*M_PI/2);
	return 1;
}
