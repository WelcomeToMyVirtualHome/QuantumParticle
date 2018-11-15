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
		const T kappa = 0;
		const T tau = 1;
		const T omega = 1;
		
		const size_t N = 100; // number of points
		const T delta_x = 1./N; 
		const size_t s_out = 10; // out.dat: info saved every 100 iterations
		
		std::vector<T> psi_r;
		std::vector<T> psi_i;
		std::vector<T> H_r;
		std::vector<T> H_i;
		std::vector<T> x;

		const char* PARAMS_FILE = "../Data/params.dat"; // parameters
		const char* R_FILE = "../Data/r.dat"; // initial positions
		const char* RHO_FILE = "../Data/rho.dat";
		const char* OUT_FILE = "../Data/out.dat"; // format "t V E_k E_tot T P"
		const char* EVAR_FILE = "../Data/E_var.dat"; // energy variance in function of tau
		
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

		T H(std::vector<T> psi, size_t k)
		{
			if(k == 0 || k == N-1)
				return 0.;
			else
				return -0.5*(psi[k+1] + psi[k-1] - 2.* psi[k])/(Utils::f_pow<T>(delta_x,2)) + kappa*(x[k] - 0.5)*psi[k]*std::sin(omega*tau);
		}

		void CalcH_i() 
		{	
			for(size_t i = 0; i < N; i++)
				H_i[i] = H(psi_i,i);
		}

		void CalcH_r()
		{
			for(size_t i = 0; i < N; i++)
				H_r[i] = H(psi_r,i);
		}
		
		void Simulation(T delta_tau, T time)
		{
			std::ofstream outputRHO;
			std::ofstream outputOUT;
			outputRHO.open(RHO_FILE, std::ios::trunc);
			outputOUT.open(OUT_FILE, std::ios::trunc);
			size_t s_d = (size_t)(time/delta_tau);
			std::cout << s_d << " " << time << " " << delta_tau << "\n";
			for(size_t s = 0; s < s_d; s++)
			{
				for(size_t i = 0; i < N; i++)
				{
					CalcH_i	();
					psi_r[i] += H_i[i]*delta_tau/2;
					CalcH_r();
					psi_i[i] -= H_r[i]*delta_tau;
					CalcH_i();
					psi_r[i] += H_i[i]*delta_tau/2;
					// std::cout << i << " " << H(psi_r,i) << " " << H(psi_i,i) << "\n";
				}
				if(s % s_out == 0)
				{
					T sq_sum = 0, sq_x_sum = 0, h_sum = 0;
					for(size_t i = 0; i < N; i++)
					{
						T rho_i = Utils::f_pow<T>(psi_r[i],2) + Utils::f_pow<T>(psi_i[i],2);
						outputRHO << x[i] << " " << rho_i << "\n\n";// << Utils::f_pow<T>(psi_r[i],2) << " " <<  Utils::f_pow<T>(psi_i[i],2) <<"\n";
						sq_sum += rho_i;
						sq_x_sum += x[i]*Utils::f_pow<T>(psi_r[i],2) + Utils::f_pow<T>(psi_i[i],2);
						h_sum += psi_r[i]*H_r[i] + psi_i[i]*H_i[i];
					}
					outputRHO << "\n";
					T N = delta_x * sq_sum;
					T X = delta_x * sq_x_sum;
					T Epsilon = delta_x * h_sum;
					outputOUT <<  N << " " << X << " " << Epsilon << "\n";
				}
			}
			outputRHO.close();	
			outputOUT.close();
		}

		
		void FlushToFiles()
		{
			std::ofstream output;
			// output.open(PARAMS_FILE, std::ios::trunc);
			// output << n<<'\n'<<m<<'\n'<<epsilon<<'\n'<<R<<'\n'<<f<<'\n'<<L<<'\n'<<a<<'\n'<<T0<<'\n'<<tau<<'\n'<<s_0<<'\n'<<s_d<<'\n'<<s_out<<'\n'<<s_xyz;
			// output.close();

			output.open(R_FILE, std::ios::trunc);
			for(size_t i = 0; i < N; i++)
			{
				output << x[i] << " " << psi_r[i] << "\n";
			}
			output.close();

			// output.open(P_FILE, std::ios::trunc);
			// for(auto p_i : p)
			// 	output << p_i;
			// output.close();
			
			// output.open(F_FILE, std::ios::trunc);
			// for(auto F_i : F)
			// 	output << F_i;
			// output.close();
		}
	};
}

#endif
