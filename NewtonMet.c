#include"NewtonMet.h"
#include"Spline.h"
#include"s21_matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

#define NEWTON_ITERATIONS 10000000

double Fshx0(double x, double* arr) { // нахождение значения производной функции
    return arr[0] * x * x + arr[1] * x + arr[2];
}


double* findZeroPoints(double* arr) // нахождение точек, где проивзодная равна 0
{
    // колво точек, где производная равна 0, записывается в ans[0] !!!!
    double a = arr[0], b = arr[1], c = arr[2];
    if (a == 0 && b != 0) {
        double* ans = calloc(2, sizeof(double));
        ans[1] = -c / b; ans[0] = 1.0;
        return ans;
    }
    else
    {
        double D = b * b - 4 * a * c;
        if (D > 0) {
            double x1, x2;
            x1 = (sqrt(D) - b) / (2.0 * a);
            x2 = (-sqrt(D) - b) / (2.0 * a);
            double* ans = calloc(3, sizeof(double));
            ans[0] = 2.0; ans[1] = x1; ans[2] = x2;
            return ans;
        }
        else if (D == 0)
        {
            double x1 = -b / (2.0 * a);
            double* ans = calloc(2, sizeof(double));
            ans[0] = 1.0; ans[1] = x1;
            return ans;
        }
        else
        {
            double* ans = calloc(1, sizeof(double));
            ans[0] = 0.0;
            return ans;
        }
    }
}





void sort(double* arr, int size)
{
    for (int i = 0; i < size; i++)
        for (int j = size - 1; j > 0; j--)
            if (arr[j] < arr[j - 1]) 
            {
                double buf = arr[j];
                arr[j] = arr[j - 1];
                arr[j - 1] = buf;
            }
}


double* make_area(double x1, double x2, double x3, double x4)
{
    double* area = calloc(2, sizeof(double));
    if ((x3 > x2) || (x4 < x1)) return area;
    else if ((x3 == x2) && (x1 != x4)) 
    {
        area[0] = x3 - 0.1; area[1] = x3 + 0.1;
    }
    else if ((x4 == x1) && (x2 != x3))
    {
        area[0] = x4 - 0.1; area[1] = x4 + 0.1;
    }
    else 
    {
        area[0] = MAX(x1, x3);
        area[1] = MIN(x2, x4);
    }
    return area;
}

Answer* NewtonSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2) {
    double a1 = f[0]; double b1 = f[1]; double c1 = f[2]; double d1 = f[3];
    double a2 = g[0]; double b2 = g[1]; double c2 = g[2]; double d2 = g[3];

    double* size = make_area(fx1, fx2, gx1, gx2);

    double f1[] = {a1, b1, c1, d1};
    double fx0[] = {a1 - a2, b1 - b2, c1 - c2, d1 - d2}; // функция формата fx - gx = 0
    double fshx0[] = {fx0[0] * 3, fx0[1] * 2, fx0[2]}; // производная этой функции

    double** answer_points = calloc(3, sizeof(double*));
    for (int i = 0; i < 3; i++) {
        answer_points[i] = calloc(2, sizeof(double));
    }
    int count_points = 0;

    double* points = findZeroPoints(fshx0);
    double* arr = calloc(points[0] + 2, sizeof(double)); // заполняем массив с концами отрезка и точками, а потом сортируем

    for (int i = 0; i < points[0] + 2; i++) {
        if (i < 2) arr[i] = size[i];
        else arr[i] = points[i - 1];
    }

    sort(arr, points[0] + 2);

    double a, b; // концы определяющие опр отрезок
    double x1, x0;

    if (points[0] != 0) // если точки есть
    {
        for (int i = 0; i < points[0] + 1; i++)
        {
            if (arr[i] < size[0] || (arr[i] == size[0] && arr[i + 1] == size[0])) continue;
            if (arr[i + 1] < size[1])
            {
                a = arr[i]; b = arr[i + 1];
            }
            else if (arr[i + 1] == size[1])
            {
                a = arr[i]; b = size[1];
            }
            else continue;

            if (Fx0(a, fx0) * Fx0(b, fx0) <= 0) // проверяем, пересекла ли функция позицию 0 в заданном отрезке
            {
                x1 = (a + b) / 2;
                x0 = x1 - 5;

                // реализация самого метода Ньютона
                while (ABS(x0 - x1) > EPS)
                {
                    x0 = x1;
                    x1 = x0 - (Fx0(x0, fx0) / Fshx0(x0, fshx0));
                }
                answer_points[count_points][0] = x1;
                answer_points[count_points][1] = Fx0(x1, f1); 
                count_points++;
                
            }
        }
    }
    if ((points[0] == 0) && (Fx0(size[0], fx0) * Fx0(size[1], fx0) <= 0))  // елси в нашем отрезке нет точек, где производная равна 0, берем весь отрезок
    {

        x1 = (size[0] + size[1]) / 2;
        x0 = x1 - 5;

        // реализация самого метода Ньютона
        while (ABS(x0 - x1) > EPS)
        {
            x0 = x1;
            x1 = x0 - (Fx0(x0, fx0) / Fshx0(x0, fshx0));
        }
        answer_points[count_points][0] = x1;
        answer_points[count_points][1] = Fx0(x1, f1); 
        count_points++;
    }

    Answer* res = (Answer*)calloc(1, sizeof(Answer));

    if (count_points == 0) {
        double dis_end = (size[0] + size[1]) / 2;
        double dis_start = dis_end - 1; 

        while (ABS(dis_end - dis_start) > EPS)
        {
            dis_start = dis_end;
            dis_end = dis_start - (((fshx0[0] * dis_start * dis_start) + (fshx0[1] * dis_start) + fshx0[2]) / ((fshx0[0] * 2 * dis_start) + fshx0[1]));
        }
        res->type = DISTANCE;
        res->distance = dis_end;      
    }
    else {
        res->type = POINT;
        res->n = count_points;
        res->point = answer_points;
    }

    free(size);
    free(points);
    free(arr);

    return res;
}

//Exit condition for newton
int condition(matrix_t* x, matrix_t* x_next, matrix_t* gradient) {
    int vars_condition = (x->matrix[0][0] - x_next->matrix[0][0]) * (x->matrix[0][0] - x_next->matrix[0][0]) + (x->matrix[0][1] - x_next->matrix[0][1]) * (x->matrix[0][1] - x_next->matrix[0][1]) > EPS;
    int grad_condition = (gradient->matrix[0][0] - gradient->matrix[0][1]) * (gradient->matrix[0][0] - gradient->matrix[0][1]) > EPS;
    return vars_condition || grad_condition;
}

//Newton method using gesse and gradient matrix
double NewtonOptimise(Spline* sp1, Spline* sp2) {
    
    int index1 = 0;
    int index2 = 0;
    int i = 0;

    matrix_t gradient;
    matrix_t gess;
    matrix_t gess_inv;
    matrix_t d;
    matrix_t x;
    matrix_t x_next;

    s21_create_matrix(2, 1, &gradient);
    s21_create_matrix(2, 2, &gess);
    s21_create_matrix(2, 2, &gess_inv);
    s21_create_matrix(2, 1, &d);
    s21_create_matrix(2, 1, &x);
    s21_create_matrix(2, 1, &x_next);

    //start points
    x.matrix[0][0] = sp1->points[0][0];
    x.matrix[0][1] = sp2->points[0][0];

    //Caluclating gradient and gesse matrix
    gradient.matrix[0][0] = FirstArgPartialDerFunction(sp1->functions[0], sp2->functions[0], x.matrix[0][0], x.matrix[0][1]);
    gradient.matrix[0][1] = SecondArgPartialDerFunction(sp1->functions[0], sp2->functions[0], x.matrix[0][0], x.matrix[0][1]);
    gess.matrix[0][0] = FirstArgPartialCommonDerivative(sp1->functions[0], sp2->functions[0], x.matrix[0][0], x.matrix[0][1]);
    gess.matrix[0][1] = AllArgReversedPartialCommonDerivative(sp1->functions[0], sp2->functions[0], x.matrix[0][0], x.matrix[0][1]);
    gess.matrix[1][0] = AllArgPartialCommonDerivative(sp1->functions[0], sp2->functions[0], x.matrix[0][0], x.matrix[0][1]);
    gess.matrix[1][1] = SecondArgPartialCommonDerivative(sp1->functions[0], sp2->functions[0], x.matrix[0][0], x.matrix[0][1]);

    //using formula x_next = x - gesse^-1 * gradient
    s21_inverse_matrix(&gess, &gess_inv);
    s21_mult_matrix(&gess_inv, &gradient, &d);
    s21_sub_matrix(&x, &d, &x_next);

    //iterations
    while ( condition(&x, &x_next, &gradient) && ++i < NEWTON_ITERATIONS) {
        x.matrix[0][0] = x_next.matrix[0][0];
        x.matrix[0][1] = x_next.matrix[0][1];

        //checking area of definition, switching functions
        for (int j = 0; j < sp1->n - 1; j++) {
            if (x.matrix[0][0] < sp1->points[j][0] && j == 0) {
                x.matrix[0][0] = (sp1->points[j][0] - x.matrix[0][0]) + sp1->points[j][0]; //mirror
                index1 = 0;
            }
            if (x.matrix[0][0] > sp1->points[j + 1][0] && j + 1 == sp1->n - 1) {
                x.matrix[0][0] = sp1->points[j][0] + x.matrix[0][0] - sp1->points[j + 1][0]; //mirror
                index1 = j;
            }
            if (x.matrix[0][0] >= sp1->points[j][0] && x.matrix[0][0] <= sp1->points[j + 1][0]) {
                index1 = j;
                break;
            }
        }
        for (int j = 0; j < sp2->n - 1; j++) {
            if (x.matrix[0][1] < sp2->points[j][0] && j == 0) {
                x.matrix[0][1] = (sp2->points[j][0] - x.matrix[0][1]) + sp2->points[j][0]; // mirror
                index2 = 0;
            }
            if (x.matrix[0][0] > sp2->points[j + 1][0] && j + 1 == sp2->n - 1) {
                x.matrix[0][1] = sp2->points[j][0] + x.matrix[0][1] - sp2->points[j + 1][0]; //mirror
                index2 = j;
            }
            if (x.matrix[0][1] >= sp2->points[j][0] && x.matrix[0][0] <= sp2->points[j + 1][0]) {
                index1 = j;
                break;
            }
        }
        //main calculations

        gradient.matrix[0][0] = FirstArgPartialDerFunction(sp1->functions[index1], sp2->functions[index2], x.matrix[0][0], x.matrix[0][1]);
        gradient.matrix[0][1] = SecondArgPartialDerFunction(sp1->functions[index1], sp2->functions[index2], x.matrix[0][0], x.matrix[0][1]);
        gess.matrix[0][0] = FirstArgPartialCommonDerivative(sp1->functions[index1], sp2->functions[index2], x.matrix[0][0], x.matrix[0][1]);
        gess.matrix[0][1] = AllArgReversedPartialCommonDerivative(sp1->functions[index1], sp2->functions[index2], x.matrix[0][0], x.matrix[0][1]);
        gess.matrix[1][0] = AllArgPartialCommonDerivative(sp1->functions[index1], sp2->functions[index2], x.matrix[0][0], x.matrix[0][1]);
        gess.matrix[1][1] = SecondArgPartialCommonDerivative(sp1->functions[index1], sp2->functions[index2], x.matrix[0][0], x.matrix[0][1]);

        s21_inverse_matrix(&gess, &gess_inv);
        s21_mult_matrix(&gess_inv, &gradient, &d);
        s21_sub_matrix(&x, &d, &x_next);

    }
    /*printf("MIN %lf %lf Dist: %lf", x_next.matrix[0][0], x_next.matrix[0][0], dist_sec_degree(sp1->functions[index1], sp2->functions[index2], x_next.matrix[0][0], x_next.matrix[0][1]));*/


    return sqrt(dist_sec_degree(sp1->functions[index1], sp2->functions[index2], x_next.matrix[0][0], x_next.matrix[0][1]));
}

//not in use
double NewtonOptimise_old(double* f, double* g, double fx1, double fx2, double gx1, double gx2, double fixed_x, int argnum) {
    double x0, x1;
    double (*func)(double* f, double* g, double x, double c);
    double (*derivative)(double* f, double* g, double x, double c);
    double l, r;
    int i = 0; // iterations
    //argnum: which argument should be used to minimise function
    if (argnum == 1) { // x1
        x0 = (fx1 + fx2) / 2;
        func = FirstArgFunction;
        derivative = FirstArgPartialDerFunction;
        l = fx1;
        r = fx2;
    }
    else { //x2
        x0 = (gx1 + gx2) / 2;
        func = SecondArgFunction;
        derivative = SecondArgPartialDerFunction;
        l = gx1;
        r = gx2;
    }

    x1 = x0 - func(f, g, x0, fixed_x) / derivative(f, g, x0, fixed_x);
    double min = func(f, g, x0, fixed_x);

    while (ABS(x0 - x1) > EPS && i < 10) {
        if (x1 < l || x1 > r) {
            x1 = (x0 + x1) / 2;
            continue;
        }
        x0 = x1;
        x1 = x0 - func(f, g, x0, fixed_x) / derivative(f, g, x0, fixed_x);
        if (min > func(f, g, x0, fixed_x))
            min = func(f, g, x0, fixed_x);
        printf("%lf %lf %lf\n", x0, func(f, g, x0, fixed_x), derivative(f, g, x0, fixed_x));

        i++;
    }
    printf("MIN %lf\n", min);
    return x1;
}

//not in use
double CoordDescent_old(double* f, double fx1, double fx2, double* g, double gx1, double gx2) {
    double x1 = (fx1 + fx2) / 2;
    double x2 = (gx1 + gx2) / 2;

    for (int i = 0; i < 1; i++) {
        if (i % 100 == 0)
            printf("%d", i);
        //fix x1, find x2
        x2 = NewtonOptimise_old(f, g, fx1, fx2, gx1, gx2, x1, 2);
        /*printf("%lf %lf\n", x1, x2);*/
        //fix x2, find x1
        x1 = NewtonOptimise_old(f, g, fx1, fx2, gx1, gx2, x2, 1);
        /*printf("%lf %lf\n", x1, x2);*/
        if (ABS(FirstArgPartialDerFunction(f, g, x1, x2)) < EPS) { //min is found
            break;
        }
    }
    return FirstArgFunction(f, g, x1, x2) * FirstArgFunction(f, g, x1, x2);
}
