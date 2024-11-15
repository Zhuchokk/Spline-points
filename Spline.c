#include"Spline.h"
#include<stdlib.h>
#include<stdio.h>

#ifndef __STDC__
#pragma warning(disable:4996)
#endif 


//Reads current file stream and returnes nearest double number, do not forget to add EOF to the file
double parse_number(FILE* fp, char end_symbol) {
	char ch;
	char number[10];
	int i = 0;

	while ((ch = fgetc(fp)) != end_symbol) {
		number[i] = ch;
		i++;
		number[i] = '\0';
	}
	return atof(number);
}

//Creates quadratic spline, using ax^2 + bx + c = y where a - given const
Spline* Constructor(char* filename, const double a) {
	FILE* fp;
	Spline* nw = (Spline*)malloc(sizeof(Spline));

	fp = fopen(filename, "r");

	nw->n = parse_number(fp, '\n');

	//Memory allocation
	nw->points = (double**)malloc(nw->n * sizeof(double*));
	for (int i = 0; i < nw->n; i++) {
		nw->points[i] = (double*)malloc(2 * sizeof(double));
	}

	nw->functions = (double**)malloc((nw->n - 1) * sizeof(double*));
	for (int i = 0; i < nw->n - 1; i++) {
		nw->functions[i] = (double*)malloc(3 * sizeof(double));
	}

	//Points parser
	for (int i = 0; i < nw->n; i++) {
		nw->points[i][0] = parse_number(fp, ' ');
		nw->points[i][1] = parse_number(fp, '\n');
	}

	//Functions creation
	for (int i = 0; i < nw->n - 1; i++) {
		double x1 = nw->points[i][0];
		double y1 = nw->points[i][1];
		double x2 = nw->points[i + 1][0];
		double y2 = nw->points[i + 1][1];

		if (x1 == x2 && y1 != y2) {
			printf("Error in creating function: x1 = x2, but y1 != y2, the contradiction in definition of a function");
			return NULL;
		}
		else if(x1 == x2 && y1 == y2){
			//the same points
			nw->functions[i][0] = a; //a
			nw->functions[i][1] = a; //b
			nw->functions[i][2] = y1 - a * x1 * x1 - nw->functions[i][1] * x1; //c
			continue;

		}

		nw->functions[i][0] = a; //a
		nw->functions[i][1] = (y2 - y1 + a * x1 * x1 - a * x2 * x2) / (x2 - x1); //b
		nw->functions[i][2] = y1 - a * x1 * x1 - nw->functions[i][1] * x1; //c
	}
	fclose(fp);
	return nw;
}