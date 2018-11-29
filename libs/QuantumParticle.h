#ifndef _QuantumParticle_h_
#define	_QuantumParticle_h_

#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <vector>
#include <limits>

#include "Utils.h"

namespace QuantumParticle
{
	template<typename T>
	class Particle
	{
	public:
		Particle()
		{ 
			psi_r.resize(N);
			psi_i.resize(N);
			x.resize(N);
			H_r.resize(N);
			H_i.resize(N);
		}
		const T L = 1;	
		const T kappa = 10;
		const std::vector<T> omegas = {3*M_PI*M_PI/2, 4*M_PI*M_PI/2, 8*M_PI*M_PI/2};
		T omega = omegas[0];
		
		const size_t N = 100; // number of points
		const T delta_x = 1./N; 
		const size_t s_out = 100; // out.dat: info saved every 10 iterations
		
		std::vector<T> psi_r;
		std::vector<T> psi_i;
		std::vector<T> H_r;
		std::vector<T> H_i;
		std::vector<T> x;

		const char* PARAMS_FILE = "../Data/params.dat"; // parameters
		const char* RHO_FILE = "../Data/rho.dat";
		const char* OUT_FILE = "../Data/out.dat"; // format "N X E"
		const char* EVAR_FILE = "../Data/e_var.dat"; // energy variance in function of tau
		const char* R_FILE = "../Data/r.dat"; // initial positions
		const char* N_FILE = "../Data/n.dat";
		const char* E_FILE = "../Data/e.dat";
		const char* X_FILE = "../Data/x.dat";
		
		void InitWaveFunction(size_t n)
		{	
			for(size_t i = 0; i < N; i++)
			{
				T x_k = i*delta_x;
				x[i] = x_k;
				psi_r[i] = std::sqrt(2)*std::sin(n*M_PI*x_k);
				psi_i[i] = 0;
			}
		}

		inline T H(std::vector<T> psi, size_t k, T tau)
		{
			if(k == 0 || k == N-1)
				return 0.;
			else
				return -0.5*(psi[k+1] + psi[k-1] - 2.* psi[k])/(Utils::f_pow<T>(delta_x,2)) + kappa*(x[k] - 0.5)*psi[k]*std::sin(omega*tau);
		}

		inline void CalculateHamiltonianR(T tau){
			for(size_t i = 0; i < N; i++)
				H_r[i] = H(psi_r,i,tau);
		}

		inline void CalculateHamiltonianI(T tau){
			for(size_t i = 0; i < N; i++)
				H_i[i] = H(psi_i,i,tau);
		}

		void Simulation(T delta_tau, T time, std::vector<size_t> nn = {1}, T n_omega = 0)
		{
			std::ofstream outputRHO, outputN, outputE, outputX, outputOUT;
			char buffer[50];
			omega = n_omega;
			for(size_t n : nn)
			{
				sprintf(buffer, "../Data/rho_n=%lu.dat", n);
				outputRHO.open(std::string(buffer), std::ios::trunc);
				outputRHO << "t rho" << "\n";
				sprintf(buffer, "../Data/n_n=%lu.dat", n);
				outputN.open(std::string(buffer), std::ios::trunc);
				sprintf(buffer, "../Data/e_n=%ld.dat", n);
				outputE.open(std::string(buffer), std::ios::trunc);
				sprintf(buffer, "../Data/x_n=%ld.dat", n);
				outputX.open(std::string(buffer), std::ios::trunc);
				sprintf(buffer, "../Data/out_n=%ld.dat", n);
				outputOUT.open(std::string(buffer), std::ios::trunc);
				size_t s_d = (size_t)(time/delta_tau);
				InitWaveFunction(n);
				T t = 0;
				for(size_t s = 0; s < s_d; s++)
				{
					CalculateHamiltonianI(t);
					for(size_t i = 0; i < N; i++)
						psi_r[i] += H_i[i]*delta_tau/2;
					t += delta_tau/2;
					CalculateHamiltonianR(t);
					for(size_t i = 0; i < N; i++)
						psi_i[i] -= H_r[i]*delta_tau;
					t += delta_tau/2;
					CalculateHamiltonianI(t);
					for(size_t i = 0; i < N; i++)
						psi_r[i] += H_i[i]*delta_tau/2;

					if(s % s_out == 0)
					{
						T sq_sum = 0, sq_x_sum = 0, h_sum = 0;
						CalculateHamiltonianR(t);
						for(size_t i = 0; i < N; i++)
						{
							T rho_i = Utils::f_pow<T>(psi_r[i],2) + Utils::f_pow<T>(psi_i[i],2);
							outputRHO << x[i] << " " << rho_i << "\n\n";
							sq_sum += rho_i;
							sq_x_sum += x[i]*Utils::f_pow<T>(psi_r[i],2) + Utils::f_pow<T>(psi_i[i],2);
							h_sum += psi_r[i]*H_r[i] + psi_i[i]*H_i[i];
						}
						outputRHO << "\n";
						T N = delta_x * sq_sum;
						T X = delta_x * sq_x_sum;
						T E = delta_x * h_sum;
						outputN << t << " " << N << "\n";
						outputE << t << " " << E << "\n";
						outputX << t << " " << X << "\n";
						outputOUT << N << " " << X << " " << E << "\n";
					}
				}
				outputRHO.close();
				outputN.close();
				outputE.close();	
				outputX.close();
				outputOUT.close();
			}
		}

		void SimulationOmegas(T delta_tau, T time, std::vector<size_t> nn = {1}, std::vector<T> omegas = {3*M_PI*M_PI/2, 4*M_PI*M_PI/2, 8*M_PI*M_PI/2})
		{
			for(size_t n : nn)
			{
				std::ofstream outputE, outputO;
				char buffer[50];
				sprintf(buffer, "../Data/omegas_n=%lu.dat",n);
				outputO.open(std::string(buffer), std::ios::trunc);
				outputO << "omega\n";	
				for(T o : omegas)
					outputO << o << "\n";	
				size_t o_index = 0;
				for(T o : omegas)
				{
					omega = o;
					sprintf(buffer, "../Data/e_n=%lu_omega=%lu.dat", n, o_index++);
					outputE.open(std::string(buffer), std::ios::trunc);
					size_t s_d = (size_t)(time/delta_tau);
					InitWaveFunction(n);
					T t = 0;
					for(size_t s = 0; s < s_d; s++)
					{
						CalculateHamiltonianI(t);
						for(size_t i = 0; i < N; i++)
							psi_r[i] += H_i[i]*delta_tau/2;
						t += delta_tau/2;
						CalculateHamiltonianR(t);
						for(size_t i = 0; i < N; i++)
							psi_i[i] -= H_r[i]*delta_tau;
						t += delta_tau/2;
						CalculateHamiltonianI(t);
						for(size_t i = 0; i < N; i++)
							psi_r[i] += H_i[i]*delta_tau/2;
						T h_sum = 0.;
						CalculateHamiltonianR(t);
						for(size_t i = 0; i < N; i++)
							h_sum += psi_r[i]*H_r[i] + psi_i[i]*H_i[i];
						outputE << t << " " << delta_x * h_sum << "\n";
					}	
					outputE.close();	
				}
				outputO.close();
			}
		}

		void SimulationResonance(T delta_tau, T time, size_t n = 1, T n_omega = 3*M_PI*M_PI/2)
		{
			std::ofstream outputO;
			char buffer[50];
			sprintf(buffer, "../Data/resonance_n=%lu.dat",n);
			outputO.open(std::string(buffer), std::ios::trunc);
			outputO << "omega_i/omega e" << "\n";
			T frac0 = 0.9;
			T frac1 = 1.1;
			size_t n_points = 10;
			T d_frac = (frac1 - frac0) / n_points;
			T initial_omega = n_omega;
			omega = frac0 * initial_omega;
			size_t s_d = (size_t)(time/delta_tau);
			for(size_t i = 0; i <= n_points; i++)
			{
				InitWaveFunction(n);
				T t = 0;
				T E_max = 0.;
				for(size_t s = 0; s < s_d; s++)
				{
					CalculateHamiltonianI(t);
					for(size_t i = 0; i < N; i++)
						psi_r[i] += H_i[i]*delta_tau/2;
					t += delta_tau/2;
					CalculateHamiltonianR(t);
					for(size_t i = 0; i < N; i++)
						psi_i[i] -= H_r[i]*delta_tau;
					t += delta_tau/2;
					CalculateHamiltonianI(t);
					for(size_t i = 0; i < N; i++)
						psi_r[i] += H_i[i]*delta_tau/2;
					T E = 0.;
					CalculateHamiltonianR(t);
					for(size_t i = 0; i < N; i++)
						E += psi_r[i]*H_r[i] + psi_i[i]*H_i[i];
					E > E_max ? E_max = E : E_max = E_max;
				}
				E_max *= delta_x;
				outputO << omega/initial_omega << " " << E_max << "\n";
				std::cout << omega/initial_omega << " " << E_max << "\n";
				frac0 += d_frac;
				omega = frac0 * initial_omega;
			}
			outputO.close();	
		}

		T SimulationEnergyVariance(T delta_tau, T time)
		{
			size_t s_d = (size_t)(time/delta_tau);
			std::vector<T> E_n(s_d);
			T E_total = 0;
			T t = 0;
			for(size_t s = 0; s < s_d; ++s)
			{	
				CalculateHamiltonianI(t);
				for(size_t i = 0; i < N; i++)
					psi_r[i] += H_i[i]*delta_tau/2;
				t += delta_tau/2;
				CalculateHamiltonianR(t);
				for(size_t i = 0; i < N; i++)
					psi_i[i] -= H_r[i]*delta_tau;
				t += delta_tau/2;
				CalculateHamiltonianI(t);
				for(size_t i = 0; i < N; i++)
					psi_r[i] += H_i[i]*delta_tau/2;
				T h_sum = 0;
				CalculateHamiltonianR(t);
				for(size_t i = 0; i < N; i++)
					h_sum += psi_r[i]*H_r[i] + psi_i[i]*H_i[i];
				E_n[s] = delta_x * h_sum;
				E_total += E_n[s];
			}	
			E_total /= s_d;
			T E_var = 0;
			for(auto e : E_n)
				E_var += Utils::f_pow<T>(E_total-e,2);
			return E_var/s_d;
		}

		void StabilityTest(size_t n, T tau1, T tau2, size_t numberOfPoints)
		{
			omega = 0;
			std::ofstream output;
			output.open(EVAR_FILE, std::ios::trunc);
			T d_tau = std::abs(tau1 - tau2)/numberOfPoints;
			T delta_tau = tau1;
			size_t N = 20000;
			for(size_t i = 0; i < numberOfPoints; ++i)
			{
				InitWaveFunction(n);
				T E = SimulationEnergyVariance(delta_tau,N*delta_tau);
				std::cout << delta_tau << " " << E << "\n";
				output << delta_tau << " " << E << "\n";
				delta_tau += d_tau;
			}
			output.close();
		}
	};
}

#endif
