#include"GradientMet.h"
#include"Spline.h"
#include<stdlib.h>

#define ITERATIONS 10000000


Answer* GradientSolve(Spline* sp1, Spline* sp2) {
	double point[2] = {sp1->points[0][0], sp2->points[0][0]};
	double point_next[2];
	double gradient[2] = { FirstArgPartialDerFunction(sp1->functions[0], sp2->functions[0], point[0], point[1]), SecondArgPartialDerFunction(sp1->functions[0], sp2->functions[0], point[1], point[0])};
	double step = 0.0001;
	int index1 = 0;
	int index2 = 0;
	Answer* res = (Answer*)calloc(1, sizeof(Answer));
	res->type = DISTANCE; // we assume distance as default
	res->distance = dist_sec_degree(sp1->functions[index1], sp2->functions[index2], point[0], point[1]);

	for (int i = 0; i < ITERATIONS; i++) {
		if (ABS(gradient[0]) <= EPS && ABS(gradient[1]) <= EPS && dist_sec_degree(sp1->functions[index1], sp2->functions[index2], point[0], point[1]) <= EPS) {
			/*printf("Intersection is found: %lf %lf", point[0], point[1]);*/
			//intersection
			res->n = 1;
			res->type = POINT;
			res->point = (double**)calloc(1, sizeof(double*));
			res->point[0] = (double*)calloc(2, sizeof(double));
			res->point[0][0] = point[0];
			res->point[0][1] = Fx0(point[0], sp1->functions[index1]);
			return res;
		}
		else if (ABS(gradient[0]) <= EPS && ABS(gradient[1]) <= EPS && dist_sec_degree(sp1->functions[index1], sp2->functions[index2], point[0], point[1]) > EPS) {
			//distance
			res->n = 1;
			res->distance = dist_sec_degree(sp1->functions[index1], sp2->functions[index2], point[0], point[1]);
			return res;
		}


		point_next[0] = point[0] - step * gradient[0];
		point_next[1] = point[1] - step * gradient[1];
		point[0] = point_next[0];
		point[1] = point_next[1];

		for (int j = 0; j < sp1->n - 1; j++) {
			if(point[0] < sp1->points[j][0] && j == 0) {
				/*printf("LAST MIN FOUND: %lf %lf\n", sp1->points[j][0], point[1]);*/ // too small point?
				point[0] = (sp1->points[j][0] - point[0]) + sp1->points[j][0]; //mirror
				index1 = 0;
			}
			if (point[0] > sp1->points[j + 1][0] && j + 1 == sp1->n - 1) {
				//printf("LAST MIN FOUND: %lf %lf\n", sp1->points[j + 1][0], point[1]); // too big point
				point[0] = sp1->points[j][0] + point[0] - sp1->points[j + 1][0]; //mirror
				index1 = j;
			}
			if (point[0] >= sp1->points[j][0] && point[0] <= sp1->points[j + 1][0]) {
				index1 = j;
				break;
			}
		}
		for (int j = 0; j < sp2->n - 1; j++) {
			if (point[1] < sp2->points[j][0] && j == 0) {
				//printf("LAST MIN FOUND: %lf %lf\n", point[0], sp2->points[j][0]); // to small point?
				point[1] = (sp2->points[j][0] - point[1]) + sp2->points[j][0]; // mirror
				index2 = 0;
			}
			if (point[0] > sp2->points[j + 1][0] && j + 1 == sp2->n - 1) {
				//printf("LAST MIN FOUND: %lf %lf\n", point[0], sp2->points[j + 1][0]); // to big point
				point[1] = sp2->points[j][0] + point[1] - sp2->points[j + 1][0]; //mirror
				index2 = j;
			}
			if (point[1] >= sp2->points[j][0] && point[0] <= sp2->points[j + 1][0]) {
				index1 = j;
				break;
			}
		}
		gradient[0] = FirstArgPartialDerFunction(sp1->functions[index1], sp2->functions[index2], point[0], point[1]);
		gradient[1] = SecondArgPartialDerFunction(sp1->functions[index1], sp2->functions[index2], point[1], point[0]);
		/*if(i % 100000 == 0)
			printf("New point %lf %lf Gradient %lf %lf Functions %d %d\n", point_next[0], point_next[1], gradient[0], gradient[1], index1, index2);*/
	}
	return res;
}
