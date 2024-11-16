#include"Spline.h"
#include<stdlib.h>
#include<stdio.h>

#ifndef __STDC__
#pragma warning(disable:4996)
#endif 


//Reads current file stream and returnes nearest double number, do not forget to add EOF to the file
double parse_number(FILE* fp, char* number, char end_symbol) {
	char ch;
	int i = 0;

	while ((ch = fgetc(fp)) != end_symbol) {
		number[i] = ch;
		i++;
		number[i] = '\0';
	}
	return atof(number);
}

//Creates cubic spline, using ax^3 + bx^2 + cx + d = y where a, b - given const
Spline* Constructor(char* filename, const double a, const double b) {
	FILE* fp;
	Spline* nw = (Spline*)malloc(sizeof(Spline));
	char number[10];

	fp = fopen(filename, "r");

	nw->n = parse_number(fp, number,'\n');

	//Memory allocation
	nw->points = (double**)malloc(nw->n * sizeof(double*));
	for (int i = 0; i < nw->n; i++) {
		nw->points[i] = (double*)malloc(2 * sizeof(double));
	}

	nw->functions = (double**)malloc((nw->n - 1) * sizeof(double*));
	for (int i = 0; i < nw->n - 1; i++) {
		nw->functions[i] = (double*)malloc(4 * sizeof(double));
	}

	//Points parser
	for (int i = 0; i < nw->n; i++) {
		nw->points[i][0] = parse_number(fp,number, ' ');
		nw->points[i][1] = parse_number(fp, number,'\n');
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
			nw->functions[i][1] = b; //b
			nw->functions[i][2] = b; //c, every is suitable
			nw->functions[i][3] = nw->functions[i][3] = y1 - a * x1 * x1 * x1 - b * x1 * x1 - nw->functions[i][2] * x1; //d
			continue;

		}

		nw->functions[i][0] = a; //a
		nw->functions[i][1] = b; //b
		nw->functions[i][2] = (y2 - a * x2 * x2 * x2 - b * x2 * x2 - y1 + a * x1 * x1 * x1 + b * x1 * x1) / (x2 - x1); //c
		nw->functions[i][3] = y1 - a * x1 * x1 * x1 - b * x1 * x1 - nw->functions[i][2] * x1; //d
	}
	fclose(fp);
	return nw;
}