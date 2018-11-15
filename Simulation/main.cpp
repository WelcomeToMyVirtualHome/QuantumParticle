#include <iostream>
#include "../libs/QuantumParticle.h"

int main(int argc, char **argv)
{
	QuantumParticle::Particle<double> p;
	p.InitWaveFunction(4);
	p.FlushToFiles();
	p.Simulation(0.0001,1);
	return 1;
}
