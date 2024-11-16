#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

#define EPS 0.00001

typedef struct Spline {
	int n;
	double** points;
	double** functions;
} Spline;

typedef struct Answer {
	int type; //0 or 1; Point = 1, Distance = 0
	int n; // number of points
	double** point;
	double distance; // min distance
} Answer;

Spline* Constructor(char* filename, const double a, const double b);

#endif