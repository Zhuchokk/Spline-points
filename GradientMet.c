#include"GradientMet.h"
#include"Spline.h"
#include<stdlib.h>

#define ITERATIONS 10000000

long double func(long double f[], long double x) { //cubic function without a module
	long double ans;
	ans = f[0] * x * x * x + f[1] * x * x + f[2] * x + f[3];
	return ans;
}

double derivative(double a, double b, double c, double x) { //derivative without a module
	double der;
	der = 3 * a * x * x + 2 * b * x + c;
	return der;
}

long double partial_derivative(long double coeff_the_fisrt_unknown[], long double x, long double the_second_unknown, long double coeff_the_second_unknown[]) { //partial derivative without a module
	long double part_der;
	part_der =
		6 * coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[0] * x * x * x * x * x +
		10 * coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[1] * x * x * x * x +
		4 * (coeff_the_fisrt_unknown[1] * coeff_the_fisrt_unknown[1] + 2 * coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[2]) * x * x * x +
		6 * (coeff_the_fisrt_unknown[0] * coeff_the_fisrt_unknown[3] + coeff_the_fisrt_unknown[1] * coeff_the_fisrt_unknown[2] - coeff_the_fisrt_unknown[0] * func(coeff_the_second_unknown, the_second_unknown)) * x * x
		+ 2 * (1 + coeff_the_fisrt_unknown[1] * coeff_the_fisrt_unknown[1] + 2 * coeff_the_fisrt_unknown[1] * coeff_the_fisrt_unknown[3] - 2 * func(coeff_the_second_unknown, the_second_unknown) * coeff_the_fisrt_unknown[1]) * x;
	+2 * (-the_second_unknown + coeff_the_fisrt_unknown[2] * coeff_the_fisrt_unknown[3] - func(coeff_the_second_unknown, the_second_unknown) * coeff_the_fisrt_unknown[2]);
	return part_der;
}

long double dist_sec_degree(long double f[], long double g[], long double x1, long double x2) { //square of the distance between the points of two splines
	long double s = (x1 - x2) * (x1 - x2) + (func(f, x1) - func(g, x2)) * (func(f, x1) - func(g, x2));
	return s;
}

Answer* GradientSolve(long double* f, long double fx1, long double fx2, long double* g, long double gx1, long double gx2) {
	//f, g - arrays of coefficients
	long double x;
	int flag = 0;
	long double distance, distance2 = -1; //distance and distance^2

	long double xf = (fx1 + fx2) / 2;
	long double xg = (gx1 + gx2) / 2;
	for (long long i = 0; i < ITERATIONS; i++) {
		if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
			x = xf;
			distance = 0;
			flag = 3;
			break;
		}
		else {
			long double xf_old = xf, xg_old = xg;
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
		for (long long i = 0; i < ITERATIONS; i++) {
			if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
				x = xf;
				distance = 0;
				flag = 3; //there's a root
				break;
			}
			else {
				long double xg_old = xg;
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
		for (long long i = 0; i < ITERATIONS; i++) {
			if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
				x = xf;
				distance = 0;
				flag = 3;
				break;
			}
			else {
				long double xf_old = xf;
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
		long double** mass1 = (long double**)calloc(1, sizeof(long double*));
		mass1[0] = (long double*)calloc(2, sizeof(long double));
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

Answer* GradientSolve_modify(long double* f, long double fx1, long double fx2, long double* g, long double gx1, long double gx2) {
	//f, g - arrays of coefficients
	long double x;
	long double distance, cheat_distance = -1; //main and auxiliary variables

	long double xf = (fx1 + fx2) / 2;
	long double xg = (gx1 + gx2) / 2;

	long double cheat_distance_min = sqrt(dist_sec_degree(f, g, xf, xg));
	long double cheat_xf = xf, cheat_xg = xg;

	for (long long i = 0; i < ITERATIONS; i++) {
		if ((xf == xg) && (func(f, xf) == func(g, xg))) { //there is a root/no root in the previous iteration or at the beginning
			x = xf;
			distance = 0;
			break;
		}
		else {
			cheat_distance = sqrt(dist_sec_degree(f, g, xf, xg));
			if (cheat_distance < cheat_distance_min) {
				cheat_distance_min = cheat_distance;
				cheat_xf = xf;
				cheat_xg = xg;
			}
			if (cheat_distance == 0) {
				x = xf;
				distance = 0;
				break;
			}
			long double xf_old = xf, xg_old = xg;
			xf = xf - EPS * partial_derivative(f, xf, xg, g);
			xg = xg - EPS * partial_derivative(g, xg, xf_old, f);

			if (xf < fx1) { //if we exit the interval: roll back to the previous step
				xf = (fx1 + xf_old) / 2;
				xg = xg_old;
				//printf("xf1: %Lf\n", xf);
			}
			if (xf > fx2) {
				xf = (fx2 + xf_old) / 2;
				xg = xg_old;
				//printf("xf2: %Lf\n", xf);
			}
			if (xg < gx1) {
				xg = (gx1 + xg_old) / 2;
				xf = xf_old;
				//printf("xg1: %Lf\n", xg);
			}
			if (xg > gx2) {
				xg = (gx2 + xg_old) / 2;
				xf = xf_old;
				//printf("xg2: %Lf\n", xg);
			}
		}
}
	distance = sqrt(dist_sec_degree(f, g, xf, xg));
	if (cheat_distance < cheat_distance_min) {
		cheat_distance_min = cheat_distance;
		cheat_xf = xf;
		cheat_xg = xg;
	}

	printf("xf: %Lf\n xg: %Lf\n x: %Lf\n distance: %Lf\n cheat_distance_min: %Lf\n cheat_xf: %Lf\n cheat_xg: %Lf\n", xf, xg, x, distance, cheat_distance_min, cheat_xf, cheat_xg);
	Answer* res = (Answer*)calloc(1, sizeof(Answer));
	if (distance == 0) {
		long double** mass1 = (long double**)calloc(1, sizeof(long double*));
		mass1[0] = (long double*)calloc(2, sizeof(long double));
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

		/*long double f_test[4] = { 1, 2, -13, 10 };
	long double g_test[4] = { 1, 2, -6.25, 1 };
	long double x_test1 = 1, x_test2 = 2, x_test3 = 0, x_test4 = 1.5;*/

		long double f_test[4] = { 1, 2, -13, 10 };
		long double g_test[4] = { 1, 2, -13, 1 };
		long double x_test1 = 0, x_test2 = 15, x_test3 = 0, x_test4 = 15;

		Answer* GradientSolve(f_test, x_test1, x_test2, g_test, x_test3, x_test4);
		Answer* GradientSolve_modify(f_test, x_test1, x_test2, g_test, x_test3, x_test4);
	}
#endif