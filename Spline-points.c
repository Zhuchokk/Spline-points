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

#ifdef MAINFILE_DEBUG


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

	Spline *sp1 = Constructor("Spline1.txt", 1, 2); // Add there variable of file name
	Spline *sp2 = Constructor("Spline2.txt", 1, 2);
	
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
	}
	default:
		printf("Error in indetifying method");
		return 0;
	}

	clock_t time = clock();

	for (int i = 0; i < sp1->n; i++) {
		for (int j = 0; j < sp2->n; j++) {
			
			Answer* ans = (*method)(sp1->functions[i], sp1->points[i][0], sp1->points[i + 1][0], sp2->functions[j], sp2->points[i][0], sp2->points[i + 1][0]);
			if (ans->type == POINT) {
				IsPointFound = 1;
				for (int k = 0; k < ans->n; k++) {
					printf("The splines intersect in (%lf, %lf)\n", ans->point[k][0], ans->point[k][1]);
				}
			}
			else if (ans->type == DISTANCE) {
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
		printf("There are no spline intersections. Min distance is %lf", MinDistance);
	}

	time = time - clock();
	double timeSpent = ((double)time) / CLOCKS_PER_SEC;

	printf("Time spent: %lf", timeSpent);

	system("PAUSE");
}

#endif