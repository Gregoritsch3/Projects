#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <iostream>
#include <functional>

using namespace std;

struct rootResults
{
  int numberOfIterations;
  double rootValue;
};

double returnZero(double, double, double);
double function1(double, double, double, double);
double derivativeOfFunction1(double, double, double, double);
double findRootByNewtonRaphson(const double&, const double&,
                                    function<double(double, double, double, double)>,
                                    function<double(double, double, double, double)>, double, double, double);

#endif
