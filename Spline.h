#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

#define EPS 0.00001

typedef enum {POINT, DISTANCE} AnswerType;
typedef enum {NEWTON, GRADIENT, QR} MethodType;

typedef struct Spline {
	int n;
	double** points;
	double** functions;
} Spline;

typedef struct Answer {
	AnswerType type; //POINT or DISTANCE
	int n; // number of points
	double** point;
	double distance; // min distance
} Answer;

Spline* Constructor(char* filename, const double a, const double b);

#endif