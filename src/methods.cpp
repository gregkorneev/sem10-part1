#include <iostream>
#include <cmath>
#include "methods.h"

using namespace std;

void printVector(const double x[], int n) {
    for (int i = 0; i < n; ++i) {
        cout << "x" << (i + 1) << " = " << x[i] << "\n";
    }
}

double residualInfinityNorm(const double A[][N], const double x[], const double b[], int n) {
    double maxR = 0.0;

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        double r = b[i] - sum;
        if (fabs(r) > maxR) {
            maxR = fabs(r);
        }
    }

    return maxR;
}

double errorNorm2(const double x[], const double x_ref[], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - x_ref[i];
        sum += d * d;
    }
    return sqrt(sum);
}

bool gaussSolve(const double A[][N], const double b[], double x[], int n) {
    double a[N][N + 1];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = A[i][j];
        }
        a[i][n] = b[i];
    }

    // прямой ход
    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double maxVal = fabs(a[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(a[i][k]) > maxVal) {
                maxVal = fabs(a[i][k]);
                pivotRow = i;
            }
        }

        if (fabs(maxVal) < 1e-12) {
            cout << "Ошибка: матрица вырождена или близка к вырожденной.\n";
            return false;
        }

        if (pivotRow != k) {
            for (int j = 0; j <= n; ++j) {
                double temp = a[k][j];
                a[k][j] = a[pivotRow][j];
                a[pivotRow][j] = temp;
            }
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = a[i][k] / a[k][k];
            for (int j = k; j <= n; ++j) {
                a[i][j] = a[i][j] - factor * a[k][j];
            }
        }
    }

    // обратный ход
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }
        if (fabs(a[i][i]) < 1e-12) {
            cout << "Ошибка: нулевой диагональный элемент на обратном ходе.\n";
            return false;
        }
        x[i] = (a[i][n] - sum) / a[i][i];
    }

    return true;
}

int jacobiSolve(const double A[][N], const double b[], double x[], int n,
                double eps, int maxIter) {
    double xOld[N];

    for (int i = 0; i < n; ++i) {
        x[i] = 0.0;
        xOld[i] = 0.0;
    }

    int iter = 0;

    while (iter < maxIter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * xOld[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        double maxDiff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = fabs(x[i] - xOld[i]);
            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }

        iter++;

        if (maxDiff < eps) {
            break;
        }

        for (int i = 0; i < n; ++i) {
            xOld[i] = x[i];
        }
    }

    return iter;
}

int gaussSeidelSolve(const double A[][N], const double b[], double x[], int n,
                     double eps, int maxIter) {
    for (int i = 0; i < n; ++i) {
        x[i] = 0.0;
    }

    int iter = 0;

    while (iter < maxIter) {
        double maxDiff = 0.0;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }

            double newXi = (b[i] - sum) / A[i][i];
            double diff = fabs(newXi - x[i]);
            if (diff > maxDiff) {
                maxDiff = diff;
            }
            x[i] = newXi;
        }

        iter++;

        if (maxDiff < eps) {
            break;
        }
    }

    return iter;
}
