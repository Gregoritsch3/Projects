using namespace std;
#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>

struct rootResults
{
  int numberOfIterations;
  double rootValue;
};

double function1(double, double, double, double);
double derivativeOfFunction1(double, double, double, double);
double returnZero(double, double, double);
double findRootByNewtonRaphson(const double&, const double&,
                                    function<double(double, double, double, double)>,
                                    function<double(double, double, double, double)>, double, double, double);

double returnZero(double Nx_0, double omega_x, double omega_y)
{
  double x0 = 5.0, precision = 1.e-6; // Možemo slobodno promijeniti x0 na 10 ili 20
  //
  double results = findRootByNewtonRaphson(x0,precision,function1,derivativeOfFunction1, Nx_0, omega_x, omega_y);
  
  return results; 
}


/*
Funkcija f(x) čiju pozitivnu nultočku tražimo
*/
double function1(double t, double Nx_0, double omega_x, double omega_y)
{
  double fac1 = -Nx_0 * (omega_x*omega_y)/(omega_y-omega_x);
  double fac2 = exp(-omega_y*t -omega_x*t) * (omega_x*exp(omega_y*t)-omega_y*exp(omega_x*t));
  
  return fac1*fac2;
}

/*
Derivacija f'(x)
*/
double derivativeOfFunction1(double t, double Nx_0, double omega_x, double omega_y)
{
  double fac1 = Nx_0 * (omega_x*omega_y)/(omega_y-omega_x);
  double fac2 = exp(-omega_y*t -omega_x*t) * (omega_x*omega_x*exp(omega_y*t)-omega_y*omega_y*exp(omega_x*t));
  
  return fac1*fac2;
}


/*
Funkcija koja traži nultočku funkcije "functionRoots" korištenjem Newton-Raphsonove metode.
Kao argumente, funkcija uzima početno (pogođeno) rješenje za nultočku x0, traženu preciznost,
samu funkciju "functionRoots" čiju nultočku tražimo, kao i funkciju njezine derivacije
"derivativeOfFunctionRoots".
*/
double findRootByNewtonRaphson(const double& x0, const double& precision,
                                    function<double(double, double, double, double)> functionRoots,
                                    function<double(double, double, double, double)> derivativeOfFunctionRoots, double Nx_0, double omega_x, double omega_y)
{
  int maxIterations = 1000; // Definiramo maksimalan broj iteracija
  
  rootResults results;
  double c = x0;
  for (int iteration = 1; iteration <= maxIterations; ++iteration)
  {
    //cout << " Rjesenje za iteraciju " << iteration << " je " << setprecision(7) << c << endl;
    //
    if(abs(functionRoots(c, Nx_0, omega_x, omega_y )) <= precision)
    {
      results.numberOfIterations = iteration;
      results.rootValue = c;
      return results.rootValue;
    }
    else
    {
      c -= functionRoots(c, Nx_0, omega_x, omega_y )/derivativeOfFunctionRoots(c, Nx_0, omega_x, omega_y );
    }
  }
  results.numberOfIterations = maxIterations;
  results.rootValue = c;
  
  return results.rootValue;
}