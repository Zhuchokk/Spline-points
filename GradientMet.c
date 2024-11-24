#include"GradientMet.h"
#include"Spline.h"

#define ITERATIONS 10000000

double func(double a, double b, double c, double d, double x) { //cubic function without a module
	double ans;
	ans = a * x * x * x + b * x * x + c * x + d;
	return ans;
}

double derivative(double a, double b, double c, double x) { //derivative without a module
	double der;
	der = 3 * a * x * x + 2 * b * x + c;
	return der;
}

double ABS(double a) { //the number module
	if (a < 0) {
		a = -a;
	}
	return a;
}

Answer* GradientSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2) {
	//f, g - arrays of coefficients
	double a_res, b_res, c_res, d_res, x1, x2, x;
	int flag;
	
	double distance = -1;
	a_res = f[0] - g[0];
	b_res = f[1] - g[1];
	c_res = f[2] - g[2];
	d_res = f[3] - g[3];
	if (fx1 < gx1) { //find the intersection of the intervals
		x1 = gx1;
	}
	else {
		x1 = fx1;
	}
	if (fx2 < gx2) {
		x2 = fx2;
	}
	else {
		x2 = gx2;
	}
	if (x2 < x1) {
		distance = -1;
	}
	else {
		if (x2 == x1) { //if the intervals have one common point
			distance = func(a_res, b_res, c_res, d_res, x2);
			ABS(distance);
		}
		else {
			x = (x1 + x2) / 2;
			double begin = func(a_res, b_res, c_res, d_res, x);
			if (begin != 0) { //if the root was found immediately
				if (begin < 0) {
					flag = -1;
				}
				else {
					flag = 1;
				}
				for (int j = 0; j < ITERATIONS; j++) { //implementation of the gradient descent method
					double x_old = x;
					x = x - EPS * flag * derivative(a_res, b_res, c_res, x);
					if ((x <= x1) || (x >= x2)) { //if the edge of the interval is reached, write the value of the function to the appropriate variable
						double help1, help2;
						help1 = func(a_res, b_res, c_res, d_res, x1);
						help2 = func(a_res, b_res, c_res, d_res, x2);
						ABS(help1);
						ABS(help2);

						if (help1 > help2) {
							distance = help2;
						}
						else {
							distance = help1;
						}
						break;
					}
					else {
						double help3 = func(a_res, b_res, c_res, d_res, x); //checking the value of the function
						if (help3 == 0) {
							distance = x;
							break;
						}
						else { //if the function changes the sign at this step, the flag value changes
							double help_old = func(a_res, b_res, c_res, d_res, x_old);
							if (help_old * help3 < 0) {
								flag = -flag;
							}
						}

					}
				}
				if (distance == -1) { //if the edge of the interval has not been reached, and the function is not zero anywhere
					distance = ABS(func(a_res, b_res, c_res, d_res, x));
				}
			}
		}
	}
	double** mass1 = (double**)calloc(1, sizeof(double*));
	double mass2[2];
	
	
	mass1[0][0] = distance;
	double** pp = &mass1;
	//printf("%lf\n %lf\n %lf\n %lf\n", x1, x2, x, distance); //not necessary for the function, but useful for tests
	Answer* res = (Answer*)calloc(1, sizeof(Answer));
	if (distance == 0) {
		res->point = pp;
	}
	else {
		res->distance = distance;
	}
	res->n = 1;

	return res;
}

#define GRADIENTMET 0
#ifndef GRADIENTMET
	int main() {
		double f_test[4] = { 2, -8, 4, 8 };
		double g_test[4] = { 1, 1, 1, 1 };
		double x_test1 = 0, x_test2 = 5, x_test3 = -3, x_test4 = 2;
		Answer* GradientSolve(f_test, x_test1, x_test2, g_test, x_test3, x_test4);
	}
#endif