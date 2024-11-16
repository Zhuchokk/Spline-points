#include<stdio.h>
#include<stdlib.h>
#include"Spline.h"
#include"NewtonMet.h"
#include"GradientMet.h"
#include"QrMet.h"

int main()
{
	Spline *sp1 = Constructor("Spline1.txt", 1, 2);
	Spline *sp2 = Constructor("Spline2.txt", 1, 2);

	for (int i = 0; i < sp1->n - 1; i++) {
		for (int g = 0; g < 4; g++) {
			printf("%lf ", sp1->functions[i][g]);
		}
		printf("\n");
	}
	
}
