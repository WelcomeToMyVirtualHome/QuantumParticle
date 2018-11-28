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
		const T kappa = 5;
		const T omega = 3*Utils::f_pow<T>(M_PI,2)/2;
		
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

		inline void Normalize()
		{
			T norm = 0;
			for(size_t i = 0; i < N; i++)
			{
				norm += Utils::f_pow<T>(psi_r[i],2) + Utils::f_pow<T>(psi_i[i],2);
			}
			for(size_t i = 0; i < N; i++)
			{
				psi_r[i] /= norm;
				psi_i[i] /= norm;
			}
		}

		inline void CalculateHamiltonianR(T tau){
			for(size_t i = 0; i < N; i++)
				H_r[i] = H(psi_r,i,tau);
		}

		inline void CalculateHamiltonianI(T tau){
			for(size_t i = 0; i < N; i++)
				H_i[i] = H(psi_i,i,tau);
		}

		void Simulation(T delta_tau, T time, std::vector<size_t> nn = {1})
		{
			std::ofstream outputRHO, outputN, outputE, outputX, outputOUT;
			char buffer[50];
			for(size_t n : nn)
			{
				sprintf(buffer, "../Data/rho_n=%lu.dat", n);
				outputRHO.open(std::string(buffer), std::ios::trunc);
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
				{
					h_sum += psi_r[i]*H_r[i] + psi_i[i]*H_i[i];
				}
				E_n[s] = delta_x * h_sum;
				E_total += E_n[s];
			}
			E_total /= s_d;
			T E_var = 0;
			for(auto e : E_n)
				E_var += Utils::f_pow<T>(E_total-e,2);
			return E_var/s_d;
		}

		void StabilityTest(size_t n, T tau1, T tau2, size_t numberOfPoints, T time)
		{
			std::ofstream output;
			output.open(EVAR_FILE, std::ios::trunc);
			T d_tau = std::abs(tau1 - tau2)/numberOfPoints;
			T delta_tau = tau1;
			for(size_t i = 0; i < numberOfPoints; ++i)
			{
				InitWaveFunction(n);
				std::cout << delta_tau << "\n";
				output << delta_tau << " " << SimulationEnergyVariance(delta_tau,time) << "\n";
				delta_tau += d_tau;
			}
			output.close();
		}

		void FlushToFiles()
		{
			std::ofstream output;
			output.open(R_FILE, std::ios::trunc);
			for(size_t i = 0; i < N; i++)
			{
				output << x[i] << " " << psi_r[i] << "\n";
			}
			output.close();
		}
	};
}

#endif
