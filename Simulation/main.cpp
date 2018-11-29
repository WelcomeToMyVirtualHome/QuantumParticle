#include <iostream>
#include "../libs/QuantumParticle.h"

int main(int argc, char **argv)
{
	QuantumParticle::Particle<double> p;
	p.Simulation(0.00002,100,{1,2,3,4},8*M_PI*M_PI/2);
	p.StabilityTest(1, std::pow(10,-8), std::pow(10,-4) * 1.2, 20);
	p.SimulationResonance(0.00002,20,1,3*M_PI*M_PI/2);
	return 1;
}
