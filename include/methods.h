#ifndef METHODS_H
#define METHODS_H

const int N = 3; // размер системы

// Вывод вектора
void printVector(const double x[], int n);

// Норма невязки (max |r_i|)
double residualInfinityNorm(const double A[][N], const double x[], const double b[], int n);

// Евклидова норма погрешности ||x - x_ref||_2
double errorNorm2(const double x[], const double x_ref[], int n);

// Метод Гаусса (прямой)
bool gaussSolve(const double A[][N], const double b[], double x[], int n);

// Метод Якоби
int jacobiSolve(const double A[][N], const double b[], double x[], int n,
                double eps, int maxIter);

// Метод Гаусса–Зейделя
int gaussSeidelSolve(const double A[][N], const double b[], double x[], int n,
                     double eps, int maxIter);

#endif // METHODS_H
