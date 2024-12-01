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


double mul_scalar(double* a, double* b, int len)
{
    double result = 0;
    for (int i = 0; i < len; i++)
        result += a[i] * b[i];
    return result;
}


double modul_vector(double* a, int len)
{
    double sum = 0;
    for (int i = 0; i < len; i++)
        sum += a[i] * a[i];
    return sqrtf(sum);
}


void mul_vector_on_val(double* vector, double val, double* res, int len)
{
    for (int i = 0; i < len; i++)
        res[i] = vector[i] * val;
}


void add_vectors(double* a, double* b, int flag, int len)
{
    if (flag == 1)
    {
        for (int i = 0; i < len; i++)
            a[i] += b[i];
    }
    else if (flag == 0)
    {
        for (int i = 0; i < len; i++)
            a[i] -= b[i];        
    }

}


void create_Q_matrix(matrix_t* sup, matrix_t* Q, int len)
{
    double vectors[3][3];
    double e_vectors[3][3];

    for (int i = 0; i < len; i++)
        for(int j = 0; j < len; j++)
            vectors[i][j] = sup->matrix[j][i];
    
    for (int i = 0; i < len; i++)
    {
        double buff_add[3];
        for (int k = 0; k < len; k++)
            buff_add[k] = 0;

        for (int j = 0; j < i; j++)
        {
            double val = mul_scalar(vectors[j], vectors[i], len) / mul_scalar(vectors[j], vectors[j], len);

            double buff_mul[3];
            for (int k = 0; k < len; k++)
                buff_mul[k] = 0;
        
            mul_vector_on_val(vectors[j], val, buff_mul, len);
            add_vectors(buff_add, buff_mul, 1, len);
        }

        add_vectors(vectors[i], buff_add, 0, len);

        double del = modul_vector(vectors[i], len);
        if (del != 0)
        {
            mul_vector_on_val(vectors[i], 1.0 / modul_vector(vectors[i], len), e_vectors[i], len);
        }
        else 
            for (int f = 0; f < len; f++)
                e_vectors[i][f] = 0;
    }

    for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++)
            Q->matrix[i][j] = e_vectors[j][i];
}

void create_R_matrix(matrix_t* sup, matrix_t* Q_trans, matrix_t* R, int len)
{
    s21_mult_matrix(Q_trans, sup, R);
    for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++)
            if (i > j) R->matrix[i][j] = 0;
}


void get_parameters(matrix_t* R, double* arr, int len)
{
    for (int i = 0; i < len; i++)
        for (int j = i; j < i + 1; j++)
            arr[i] = R->matrix[i][j];
}


double AbS(double val)
{
    if (val >= 0) return val;
    return -val;
}


double max_sub_abs(double* a, double* b, int len)
{
    double max = -999999999;
    for (int i = 0; i < len; i++)
        if (AbS(a[i] - b[i]) > max) max = AbS(a[i] - b[i]);
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

    Answer* res = (Answer*)calloc(1, sizeof(Answer));

    int count_answers = 0;
    int len = 0;

    if ((function[0] != 0) && (function[0] != 1))
    {
        double del = function[0];
        for (int i = 0; i < 4; i++)
            function[i] /= del;
    }
    else if ((function[1] != 0) && (function[1] != 1) && (function[0] == 0))
    {
        double del = function[1];
        for (int i = 1; i < 4; i++)
            function[i] /= del;
    }
    else if ((function[2] != 0) && (function[2] != 1) && (function[0] == 0) && (function[1] == 0))
    {
        if (function[3] == 0)
        {
            answer_points[0][0] = 0; answer_points[0][1] = f[3];
            if ((0 >= size[0]) && (0 <= size[1]))
            {
                res->type = POINT;
                res->n = 1;
                res->point = answer_points;
            }
            return res;
        }
        else if (function[2] == 0) return res;
        else
        {
            answer_points[0][0] = -function[3] / function[2];
            double x = answer_points[0][0];
            if ((x >= size[0]) && (x <= size[1]))
            {
                answer_points[0][1] = x * x * x * f[0] + x * x * f[1] + x * f[2] + f[3];
                res->type = POINT;
                res->n = 1;
                res->point = answer_points;
            }
            return res;
        } 
    }
    else if ((function[0] == 0) && (function[1] == 1) && (function[2] == 0) && (function[3] == 0))
    {
        return res;
    }
    if (function[0] != 0) len = 3;
    else len = 2;

    double answer_start[3];
    double answer_end[3];

    matrix_t support;
    matrix_t Q;
    matrix_t Q_trans;
    matrix_t R;

    s21_create_matrix(len, len, &support);
    s21_create_matrix(len, len, &Q);
    s21_create_matrix(len, len, &Q_trans);
    s21_create_matrix(len, len, &R);

    for (int i = 0; i < len; i++) // создаем сопровождающую матрицу
    {
        for(int j = 0; j < len; j++)
        {
            if (j == (i - 1)) support.matrix[i][j] = 1;
            else if (j == len - 1) support.matrix[i][j] = -function[3 - i];
            else support.matrix[i][j] = 0;
        }
    }
    
    create_Q_matrix(&support, &Q, len);
    s21_transpose(&Q, &Q_trans);
    create_R_matrix(&support, &Q_trans, &R, len);
    get_parameters(&R, answer_end, len);

    // задаем 100% шанс проийти условие while
    for (int i = 0; i < len; i++)
        answer_start[i] = answer_end[i] - 1;

    // QR метод
    while(max_sub_abs(answer_end, answer_start, len) > EPS)
    {
        for (int i = 0; i < len; i++)
            answer_start[i] = answer_end[i];
        s21_mult_matrix(&R, &Q, &support);
        create_Q_matrix(&support, &Q, len);
        s21_transpose(&Q, &Q_trans);
        create_R_matrix(&support, &Q_trans, &R, len);
        get_parameters(&R, answer_end, len);     
    }

    // записываем корни
    for (int i = 0; i < len; i++)
    {
        if ((answer_end[i] >= size[0]) && (answer_end[i] <= size[1]))
        {
            double x = answer_end[i];
            double y = x * x * x * (f[0] - g[0]) + x * x * (f[1] - g[1]) + x * (f[2] - g[2]) + f[3] - g[3];
            if ((y > 0.0001) || (y < -0.0001))
                answer_end[i] = -answer_end[i];
            
            x = answer_end[i];
            y = x * x * x * (f[0] - g[0]) + x * x * (f[1] - g[1]) + x * (f[2] - g[2]) + f[3] - g[3];

            if ((y <= 0.0001) && (y >= -0.0001))
            {
                answer_points[count_answers][0] = x;
                answer_points[count_answers][1] = x * x * x * f[0] + x * x * f[1] + x * f[2] + f[3];
                count_answers++;
            }
        }
    }


    if (count_answers != 0)
    {
        res->type = POINT;
        res->n = count_answers;
        res->point = answer_points;
        return res;
    }
    return res;
}