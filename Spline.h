#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

#define EPS 0.00001

typedef enum {POINT, DISTANCE} AnswerType;
typedef enum {NEWTON, GRADIENT, QR} MethodType;

typedef struct Spline {
	int n; // number of points
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

double ABS(double x);

double Fx0(double x, double* arr);

double FirstPartialDerivative(double* f, double x);
double SecondPartialDerivative(double* f, double x);

double FirstArgPartialDerFunction(double* f, double* g, double x1, double c);
double SecondArgPartialDerFunction(double* f, double* g, double x2, double c);

double FirstArgFunction(double* f, double* g, double x1, double c);
double SecondArgFunction(double* f, double* g, double x2, double c);

double FirstArgPartialCommonDerivative(double* f, double* g, double x1, double c);
double SecondArgPartialCommonDerivative(double* f, double* g, double x2, double c);

double AllArgPartialCommonDerivative(double* f, double* g, double x1, double x2);
double AllArgReversedPartialCommonDerivative(double* f, double* g, double x1, double x2);

long double dist_sec_degree(long double f[], long double g[], long double x1, long double x2);
#endif