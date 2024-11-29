#include"GradientMet.h"
#include"Spline.h"
#include<stdlib.h>

#define ITERATIONS 10000000

double func(double f[], double x) { //cubic function without a module
	double ans;
	ans = f[0] * x * x * x + f[1] * x * x + f[2] * x + f[3];
	return ans;
}

double derivative(double a, double b, double c, double x) { //derivative without a module
	double der;
	der = 3 * a * x * x + 2 * b * x + c;
	return der;
}

double partial_derivative(double coeff_the_fisrt_unknown[], double x, double the_second_unknown, double coeff_the_second_unknown[]) { //partial derivative without a module
	double part_der;
	part_der =
		6 * coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[0] * x * x * x * x * x +
		10 * coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[1] * x * x * x * x +
		4 * (coeff_the_fisrt_unknown[1] * coeff_the_fisrt_unknown[1] + 2 * coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[2]) * x * x * x +
		6 * (coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[3] + coeff_the_fisrt_unknown[1] * coeff_the_fisrt_unknown[2] - coeff_the_fisrt_unknown[0] * func(coeff_the_second_unknown, the_second_unknown)) * x * x
		+ 2 * coeff_the_fisrt_unknown[1] * (coeff_the_fisrt_unknown[1] + 2 * coeff_the_fisrt_unknown[3] - 2 * func(coeff_the_second_unknown, the_second_unknown)) * x;
	+2 * coeff_the_fisrt_unknown[2] * (coeff_the_fisrt_unknown[3] - func(coeff_the_second_unknown, the_second_unknown));
	return part_der;
}

double dist_sec_degree(double f[], double g[], double x1, double x2) { //square of the distance between the points of two splines
	double s = (x1 - x2) * (x1 - x2) + (func(f, x1) - func(g, x2)) * (func(f, x1) - func(g, x2));
	return s;
}

Answer* GradientSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2) {
	//f, g - arrays of coefficients
	double x;
	int flag = 0;
	double distance, distance2 = -1; //distance and distance^2

	double xf = (fx1 + fx2) / 2;
	double xg = (gx1 + gx2) / 2;
	for (int i = 0; i < ITERATIONS; i++) {
		if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
			x = xf;
			distance = 0;
			flag = 3;
			break;
		}
		else {
			double xf_old = xf, xg_old = xg;
			xf = xf - EPS * partial_derivative(f, xf, xg, g);
			xg = xg - EPS * partial_derivative(g, xg, xf_old, f);
			if (xf < fx1) { //if we exit the interval: roll back to the previous step and exit the loop
				xf = fx1;
				xg = xg_old;
				flag = 1;
				break;
			}
			if (xf > fx2) {
				xf = fx2;
				xg = xg_old;
				flag = -1;
				break;
			}
			if (xg < gx1) {
				xg = gx1;
				xf = xf_old;
				flag = 2;
				printf("suka %lf\n", xg_old);
				break;
			}
			if (xg > gx2) {
				xg = gx2;
				xf = xf_old;
				flag = -2;
				break;
			}
		}
	}
	if ((flag == 1) || (flag == -1)) { //fix xf and change xg
		for (int i = 0; i < ITERATIONS; i++) {
			if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
				x = xf;
				distance = 0;
				flag = 3; //there's a root
				break;
			}
			else {
				double xg_old = xg;
				xg = xg - EPS * partial_derivative(g, xg, xf, f);
				if (xg < gx1) { //when exiting and exiting the second interval, the answer is the distance between the extreme points
					xg = gx1;
					distance = sqrt(dist_sec_degree(f, g, xf, xg));
					flag = 4;
					printf("xf fixed, xg out down\n");
					break;
				}
				if (xg > gx2) {
					xg = gx2;
					distance = sqrt(dist_sec_degree(f, g, xf, xg));
					flag = 4;
					printf("xf fixed, xg out up\n");
					break;
				}
			}
		}
	}
	if ((flag == 2) || (flag == -2)) { //fix xg and change xf
		for (int i = 0; i < ITERATIONS; i++) {
			if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
				x = xf;
				distance = 0;
				flag = 3;
				break;
			}
			else {
				double xf_old = xf;
				xf = xf - EPS * partial_derivative(f, xf, xg, g);
				if (xf < fx1) { //when exiting and exiting the second interval, the answer is the distance between the extreme points
					xf = fx1;
					distance = sqrt(dist_sec_degree(f, g, xf, xg));
					flag = 4;
					printf("xg fixed, xf out down\n");
					break;
				}
				if (xf > fx2) {
					xf = fx2;
					distance = sqrt(dist_sec_degree(f, g, xf, xg));
					printf("xg fixed, xf out up\n");
					flag = 4;
					break;
				}
			}
		}
	}
	if ((flag == 0) || (flag == 1) || (flag == -1) || (flag == 2) || (flag == -2) || (flag == 4)) { //no exits from the interval/one exit from the interval
		distance = sqrt(dist_sec_degree(f, g, xf, xg));
		if (distance == 0) {
			x = xf;
		}
	}

	/*
	flag = +-1 - out of interval xf
	flag = +-2 - out of interval õg
	flag = 3 - found a root
	flag = 4 - out of both intervals
	*/
	
	//printf("%lf\n %lf\n %lf\n %lf\n flag: %d\n", xf, xg, x, distance, flag); //not necessary for the function, but useful for tests
	Answer* res = (Answer*)calloc(1, sizeof(Answer));
	if (distance == 0) {
		double** mass1 = (double**)calloc(1, sizeof(double*));
		mass1[0] = (double*)calloc(2, sizeof(double));
		mass1[0][0] = x;
		mass1[0][1] = func(f, x);

		res->type = POINT;
		res->n = 1;
		res->point = mass1;
	}
	else {
		res->type = DISTANCE;
		res->distance = distance;
	}
	

	return res;
}

#define GRADIENTMET 0
#ifndef GRADIENTMET
	int main() {
		//two options for input data for tests:

		/*double f_test[4] = { 1, 2, -13, 10 };
	double g_test[4] = { 1, 2, -6.25, 1 };
	double x_test1 = 1, x_test2 = 2, x_test3 = 0, x_test4 = 1.5;*/

		double f_test[4] = { 1, 2, -13, 10 };
		double g_test[4] = { 1, 2, -13, 1 };
		double x_test1 = 0, x_test2 = 15, x_test3 = 0, x_test4 = 15;

		Answer* GradientSolve(f_test, x_test1, x_test2, g_test, x_test3, x_test4);
	}
#endif