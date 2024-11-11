#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

struct Spline {
	int n;
	double** points;
	double** functions;
};

struct Spline* Constructor(int n, double** points, double** functions);

#endif