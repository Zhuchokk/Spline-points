#include<stdio.h>
#include<stdlib.h>
#include"Spline.h"
#include"NewtonMet.h"
#include"GradientMet.h"
#include"QrMet.h"

int main()
{
	Spline *sp1 = Constructor("Spline1.txt", 1);
	Spline *sp2 = Constructor("Spline2.txt", 1);
}
