#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"Spline.h"
#include"NewtonMet.h"
#include"GradientMet.h"
#include"QrMet.h"

#ifndef __STDC__
#pragma warning(disable:4996)
#endif 

#define MAINFILE_DEBUG 1
#if MAINFILE_DEBUG


int main()
{
	char* filename1[100];
	char* filename2[100];
	int IsPointFound = 0;
	double MinDistance = -1; // -1 means tha min distance isn't defined
	MethodType mtype;
	Answer* (*method)(double* f, double fx1, double fx2, double* g, double gx1, double gx2);

	printf("Rule for files: points must be sorted\n");
	printf("Enter the filename for the 1 Spline\n");
	scanf("%s", &filename1);

	printf("Enter the filename for the 2 Spline\n");
	scanf("%s", &filename2);

	printf("Choose the method(0 - Newton, 1 - Gradient, 2 - Accompanying matrix + QR)\n");
	scanf("%d", &mtype);

	Spline *sp1 = Constructor(filename1, 1, 2);
	Spline *sp2 = Constructor(filename2, 1, 2);
	
	//Printing points
	/*for (int i = 0; i < sp1->n - 1; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%lf", sp1->functions[i][j]);
		}
		printf("\n");
	}
	printf("\nSecond\n");
	for (int i = 0; i < sp2->n - 1; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%lf", sp2->functions[i][j]);
		}
		printf("\n");
	}*/

	switch (mtype)
	{
	case NEWTON: {
		method = &NewtonSolve;
		break;
	}
	case GRADIENT: {
		method = &GradientSolve;
		break;
	}
	case QR: {
		method = &QrSolve;
		break;
	}
	default:
		printf("Error in indetifying method");
		return 0;
	}

	clock_t time = clock();

	for (int i = 0; i < sp1->n - 1; i++) {
		for (int j = 0; j < sp2->n - 1; j++) {

			Answer* ans = (*method)(sp1->functions[i], sp1->points[i][0], sp1->points[i + 1][0] - EPS * (i != sp1->n - 2), sp2->functions[j], sp2->points[j][0], sp2->points[j + 1][0] - EPS * (j != sp2->n - 2));
			if (ans->type == POINT) {
				IsPointFound = 1;
				for (int k = 0; k < ans->n; k++) {
					printf("The splines intersect in (%lf, %lf)\n", ans->point[k][0], ans->point[k][1]);
				}
			}
			else if (ans->type == DISTANCE) {
				if (mtype == QR) { continue; } // QR can't find distance, skip
				if (MinDistance == -1) {
					MinDistance = ans->distance;
				}
				else {
					if (MinDistance > ans->distance) {
						MinDistance = ans->distance;
					}
				}
			}
			else {
				printf("Answer type error");
				return 0;
			}
		}
	}
	if (!IsPointFound) {
		double MinDistanceOptimised;
		if (mtype == NEWTON) {
			MinDistanceOptimised = NewtonOptimise(sp1, sp2);
		}
		else {
			MinDistanceOptimised = GradientOptimise(sp1, sp2);
		}
		if (MinDistanceOptimised > MinDistance &&MinDistance != -1) { // mindistance is smaller and defined
			MinDistanceOptimised = MinDistance;
		}
		printf("There are no spline intersections. Min distance is %lf\n", MinDistanceOptimised);
	}

	time = clock() - time;
	double timeSpent = ((double)time) / CLOCKS_PER_SEC;

	printf("Time spent: %lf\n", timeSpent);

	system("PAUSE");
}

#endif