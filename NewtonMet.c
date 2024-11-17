#include"NewtonMet.h"
#include"Spline.h"


Answer* NewtonSolve(double* f, double fx1, double fx2, double* g, double gx1, double gx2) {}

#define DEBAG_NEWTONMET 1
#ifdef DEBAG_NEWTONMET

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.000001


double Fx0(double x, double* arr) { // нахождение значения функции
    double ans = arr[0] * x * x * x + arr[1] * x * x + arr[2] * x + arr[3];
    return ans;
}
    

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


double ABS(double x)  // модуль числа
{
    if (x < 0) return -x;
    return x;
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

int main()
{
    double size[2] = {-3.0, 5.0}; // промежуток, где существуют одновременно два сплайна

    double fx0[] = {3.0, -6.0, -10.0, -10.0}; // функция формата fx - gx = 0
    double fshx0[] = {9.0, 12.0, -10.0}; // производная этой функции

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

            printf("%llf %llf\n", a, b);

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
                printf("x1 = %llf\n", x1);
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
        printf("x1 = %llf\n", x1);
    }


    scanf("%d");
    return 0;
}
#endif