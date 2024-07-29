#include "newtonraphson.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <random>


using namespace std;
using std::vector;


vector<vector<double>> calculateRK4_Decay(double, double , const double& , const double& , const double&, vector<double> , vector<double>, vector<double>, const double& N);
double calculate_f_Nx(const double&, const double&);
double calculate_f_Ny(const double&, const double&, const double&, const double&);
double calculate_random_number();
vector<vector<double>> calculate_MonteCarlo_Decay(const double&, const double&, const double&, const double&, const int&, const double&, const double&, const int&);


int main(){
  double Nx_0 = 10; // !!
  double Ny_0 = 0;
  double omega_x=1.0/7.2; //dana^-1
  double omega_y=1.0/200; //dana^-1
    
  double tmin=0;
  double tmax=1000;
  
  double N=10000; //broj RK4 koraka !!
  double M=10000; //broj MC koraka !!
  int cycles=1;
  
  double h=(tmax-tmin)/N;
  
  double Nx=Nx_0;
  double Ny=Ny_0;
    
  vector<double> Nx_t;
  vector<double> Ny_t;
  vector<double> t_i;
  
  Nx_t.push_back(Nx_0);
  Ny_t.push_back(Ny_0);
  t_i.push_back(0);

  vector<vector<double>> RK4_res = calculateRK4_Decay(Nx, Ny, omega_x, omega_y, h, t_i, Nx_t, Ny_t, N);
  vector<vector<double>> MC_res = calculate_MonteCarlo_Decay(Nx, Ny, tmax, tmin, cycles, omega_x, omega_y, M);
  string delim = "";
  
  /*
  //print Nx
  cout << "Nx(t)_RK4: " << endl;
  for(auto&& item : RK4_res[0])
      {
        cout << delim << item;
        delim = ", ";
      }
  cout << endl;
  
  //print Ny
  cout << "Ny(t)_RL4: " << endl;
  for(auto&& item : RK4_res[1])
      {
        cout << delim << item;
        delim = ", ";
      }
  cout << endl;
  
  
  //print t
  cout << "t: " << endl;
  for(auto&& item : res[2])
      {
        cout << delim << item;
      }
  cout << endl;
  
  
  
  //print MC Decay N_x
  cout << "Nx_MC(t): " << endl;
  
  for(auto&& item : MC_res[0])
      {
        cout << delim << item;
        delim = ", ";
      }
  cout << endl;
  cout << "-----------------------------------------------------------------------------"<< endl;
    
  //print MC Decay N_y
  cout << "Ny_MC(t): " << endl;
  
  for(auto&& item : MC_res[1])
      {
        cout << delim << item;
        delim = ", ";
      }
  cout << endl;
  cout << "-----------------------------------------------------------------------------" << endl;
  
  //print alpha_particleMC(t)
  cout << "alpha_particle_MC(t): " << endl;
  
  for(auto&& item : MC_res[2])
      {
        cout << delim << item;
        delim = ", ";
      }
  cout << endl;
  cout << "-----------------------------------------------------------------------------" << endl;
  */
  
  //Nultočke druge derivacije od N_alpha(t) u svrhu pronalaska kritičnog vremena t_crit
  double t_crit = returnZero(Nx_0, omega_x, omega_y);
  cout << t_crit << endl;
  
  return 0;
}


vector<vector<double>> calculateRK4_Decay(double Nx, double Ny, const double& omega_x, const double& omega_y, const double& h, vector<double> t_i, vector<double> Nx_t, vector<double> Ny_t, const double& N)
{
  double k1_Nx, k2_Nx, k3_Nx, k4_Nx;  
  double k1_Ny, k2_Ny, k3_Ny, k4_Ny;
  for(int i=1; i<=N; ++i)
  {
           
    k1_Nx = h*calculate_f_Nx(Nx, omega_x);
    k2_Nx = h*calculate_f_Nx(Nx+k1_Nx/2, omega_x);
    k3_Nx = h*calculate_f_Nx(Nx+k2_Nx/2, omega_x);
    k4_Nx = h*calculate_f_Nx(Nx+k3_Nx, omega_x);

    k1_Ny = h*calculate_f_Ny(Nx, Ny, omega_x, omega_y);
    k2_Ny = h*calculate_f_Ny(Nx, Ny+k1_Ny/2, omega_x, omega_y);
    k3_Ny = h*calculate_f_Ny(Nx, Ny+k2_Ny/2, omega_x, omega_y);
    k4_Ny = h*calculate_f_Ny(Nx, Ny+k3_Ny, omega_x, omega_y);
    
    Nx+=1.0/6.0*(k1_Nx+2*k2_Nx+2*k3_Nx+k4_Nx);
    Ny+=1.0/6.0*(k1_Ny+2*k2_Ny+2*k3_Ny+k4_Ny);
    
    t_i.push_back(i*h);
    Nx_t.push_back(Nx);
    Ny_t.push_back(Ny);
  }
  vector<vector<double>> results;
  results.push_back(Nx_t);
  results.push_back(Ny_t);
  results.push_back(t_i);
  
  return results;
}


double calculate_f_Nx(const double& Nx, const double& omega_x)
{
  return -omega_x*Nx;
}


double calculate_f_Ny(const double& Nx, const double& Ny, const double& omega_x, const double& omega_y)
{
  return -omega_y*Ny+omega_x*Nx;
}


double calculate_random_number()
{
  thread_local std::random_device rd;
  thread_local std::mt19937 gen(rd());
  std::uniform_real_distribution <> dist(0.0, 1.0);
  return dist(gen);
}


vector<vector<double>> calculate_MonteCarlo_Decay(const double& Nx_0, const double& Ny_0, const double& tmax, const double& tmin, const int& cycles, const double& omega_x, const double& omega_y, const int& M)
{
  vector<double> ncumulative(M+1);
  vector<double> ncumulative2(M+1);
  vector<double> nalpha(M+1);

  vector<vector<double>> ncycles;
  double h = (tmax-tmin)/M;
    
  for(int i = 1; i <= cycles; ++i)
  {
    double n_unstable = Nx_0;
    double n_unstable2 = Ny_0;
    double alpha_particle = 0;
    
    ncumulative[0]+=Nx_0/cycles;
    ncumulative2[0]+=Ny_0/cycles;
    nalpha[0]+=alpha_particle/cycles;

    for(int j = 1; j<=M; j++)
    {
      //double t = j*h;
      double decay_probability_x = omega_x*h;
      double decay_probability_y = omega_y*h;
      
      double particle_limit = n_unstable;
      double particle_limit2 = n_unstable2;

      for(int nth_particle = 1; nth_particle <= particle_limit; nth_particle++)
      {
        //random number generator
        double x = calculate_random_number();
        
        //Monte Carlo condition
        if(x <= decay_probability_x)
        {
          n_unstable=n_unstable-1;
          n_unstable2=n_unstable2+1;
        }
      }
      
      for(int nth_particle2 = 1; nth_particle2 <= particle_limit2; nth_particle2++)
      {
        //random number generator
        double y = calculate_random_number();
        
        //Monte Carlo condition
        if(y <= decay_probability_y)
        {
          n_unstable2=n_unstable2-1;
          alpha_particle+=1;
        }
      }
      ncumulative[j]+=n_unstable/cycles;
      ncumulative2[j]+=n_unstable2/cycles;
      nalpha[j]+=alpha_particle/cycles;
    }
  }
  
  vector<vector<double>> results;
  results.push_back(ncumulative);
  results.push_back(ncumulative2);
  results.push_back(nalpha);


  return results;
}
