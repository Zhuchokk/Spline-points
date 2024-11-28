#include"QrMet.h"
#include"Spline.h"
#include "s21_matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

#define EPS 0.00001


double* make_size(double x1, double x2, double x3, double x4)
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


double mul_scalar(double* a, double* b)
{
    double result = 0;
    for (int i = 0; i < 3; i++)
        result += a[i] * b[i];
    return result;
}


double modul_vector(double* a)
{
    return sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}


void mul_vector_on_val(double* vector, double val, double* res)
{
    for (int i = 0; i < 3; i++)
        res[i] = vector[i] * val;
}


void add_vectors(double* a, double* b, int flag)
{
    if (flag == 1)
    {
        for (int i = 0; i < 3; i++)
            a[i] += b[i];
    }
    else if (flag == 0)
    {
        for (int i = 0; i < 3; i++)
            a[i] -= b[i];        
    }

}


void create_Q_matrix(matrix_t* sup, matrix_t* Q)
{
    double vectors[3][3];
    double e_vectors[3][3];

    for (int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            vectors[i][j] = sup->matrix[j][i];
    
    for (int i = 0; i < 3; i++)
    {
        double buff_add[3] = {0, 0, 0};
        for (int j = 0; j < i; j++)
        {
            double val = mul_scalar(vectors[j], vectors[i]) / mul_scalar(vectors[j], vectors[j]);

            double buff_mul[3] = {0, 0, 0};
            mul_vector_on_val(vectors[j], val, buff_mul);
            add_vectors(buff_add, buff_mul, 1);
        }

        add_vectors(vectors[i], buff_add, 0);

        mul_vector_on_val(vectors[i], 1.0 / modul_vector(vectors[i]), e_vectors[i]);
    }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Q->matrix[i][j] = e_vectors[j][i];
}

void create_R_matrix(matrix_t* sup, matrix_t* Q_trans, matrix_t* R)
{
    s21_mult_matrix(Q_trans, sup, R);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (i > j) R->matrix[i][j] = 0;
}


void get_parameters(matrix_t* R, double* arr)
{
    for (int i = 0; i < 3; i++)
        for (int j = i; j < i + 1; j++)
            arr[i] = R->matrix[i][j];
}


double ABS(double val)
{
    if (val >= 0) return val;
    return -val;
}


double max_sub_abs(double* a, double* b)
{
    double max = -999999999;
    for (int i = 0; i < 3; i++)
        if (ABS(a[i] - b[i]) > max) max = ABS(a[i] - b[i]);
    return max;
}


Answer* QrSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2) 
{
    double* size = make_size(fx1, fx2, gx1, gx2);
    double function[] = {f[0] - g[0], f[1] - g[1], f[2] - g[2], f[3] - g[3]};

    double** answer_points = calloc(3, sizeof(double*));
    for (int i = 0; i < 3; i++) {
        answer_points[i] = calloc(2, sizeof(double));
    }

    int count_answers = 0;

    if (function[0] != 1.0) // делаем многочлен унитарным
        for (int i = 1; i < 4; i++)
            function[i] /= function[0]; 
 
    double answer_start[3];
    double answer_end[3];

    matrix_t support;
    matrix_t Q;
    matrix_t Q_trans;
    matrix_t R;

    s21_create_matrix(3, 3, &support);
    s21_create_matrix(3, 3, &Q);
    s21_create_matrix(3, 3, &Q_trans);
    s21_create_matrix(3, 3, &R);

    for (int i = 0; i < 3; i++) // создаем сопровождающую матрицу
    {
        for(int j = 0; j < 3; j++)
        {
            if (j == (i - 1)) support.matrix[i][j] = 1;
            else if (j == 2) support.matrix[i][j] = -function[3 - i];
            else support.matrix[i][j] = 0;
        }
    }

    create_Q_matrix(&support, &Q);
    s21_transpose(&Q, &Q_trans);
    create_R_matrix(&support, &Q_trans, &R);
    get_parameters(&R, answer_end);

    // задаем 100% шанс проийти условие while
    for (int i = 0; i < 3; i++)
        answer_start[i] = answer_end[i] - 1;

    // QR метод
    while(max_sub_abs(answer_end, answer_start) > EPS)
    {
        for (int i = 0; i < 3; i++)
            answer_start[i] = answer_end[i];
        s21_mult_matrix(&R, &Q, &support);
        create_Q_matrix(&support, &Q);
        s21_transpose(&Q, &Q_trans);
        create_R_matrix(&support, &Q_trans, &R);
        get_parameters(&R, answer_end);     
    }

    // записываем корни
    for (int i = 0; i < 3; i++)
    {
        double rev = -answer_end[i];
        // костыль, который исправляет баг с корнями (почему-то выводит отрицательные корни без минуса)
        if (rev * rev * rev * function[0] + rev * rev * function[1] + rev * function[2] + function[3] < 0.001) answer_end[i] = rev; 
        if ((answer_end[i] >= size[0]) && (answer_end[i] <= size[1]))
        {
            double x = answer_end[i] // создал для читабельности кода
            answer_points[count_answers][0] = x;
            answer_points[count_answers][1] = x * x * x * f[0] + x * x * f[1] + x * f[2] + f[3];
            count_answers++;
        }
    }
    if (count_answers == 0)
    {
        struct Answer* res = (Answer*)calloc(1, sizeof(Answer));
        res->type = POINT;
        res->n = count_answers;
        res->point = answer_points; 
    }

    free(size);
    return res;
}