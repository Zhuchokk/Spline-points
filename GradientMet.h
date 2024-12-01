#include"Spline.h"

#ifndef GRADIENT_INCLUDED
#define GRADIENT_INCLUDED

double GradientOptimise(Spline* sp1, Spline* sp2);
Answer* GradientSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2);

#endif