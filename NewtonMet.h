#include"Spline.h"

#ifndef NEWTON_INCLUDED
#define NEWTON_INCLUDED

Answer* NewtonSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2);
double NewtonOptimise(Spline* sp1, Spline* sp2);

#endif