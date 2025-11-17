#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>

#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #include <sys/types.h>
    #define MKDIR(path) mkdir(path, 0755)
#endif

#include "../include/benchmark.h"

using namespace std;

// Простая обёртка для создания папки csv
static void ensureCsvDirBenchmark() {
    int result = MKDIR("csv");
    (void)result; // игнорируем код возврата (если папка уже есть)
}

// Индексация элемента A[i][j] в плоском массиве
static inline int idx(int i, int j, int n) {
    return i * n + j;
}

// Генерация диагонально преобладающей матрицы A (n x n) и вектора b (размер n)
static void generateDiagonalDominantSystem(double* A, double* b, int n) {
    for (int i = 0; i < n; ++i) {
        double rowSum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double val = ((double)rand() / RAND_MAX) * 0.2 - 0.1; // [-0.1; 0.1]
            A[idx(i, j, n)] = val;
            rowSum += fabs(val);
        }
        A[idx(i, i, n)] = rowSum + 1.0; // диагональное преобладание
        b[i] = 1.0;                      // правая часть просто 1
    }
}

// --- Динамический метод Гаусса (без выбора главного элемента) ---
static void gaussSolveDyn(double* A, double* b, double* x, int n) {
    // прямой ход
    for (int k = 0; k < n; ++k) {
        double pivot = A[idx(k, k, n)];
        if (fabs(pivot) < 1e-12) continue; // грубая защита

        for (int i = k + 1; i < n; ++i) {
            double factor = A[idx(i, k, n)] / pivot;
            for (int j = k; j < n; ++j) {
                A[idx(i, j, n)] -= factor * A[idx(k, j, n)];
            }
            b[i] -= factor * b[k];
        }
    }

    // обратный ход
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[idx(i, j, n)] * x[j];
        }
        double diag = A[idx(i, i, n)];
        if (fabs(diag) < 1e-12) {
            x[i] = 0.0;
        } else {
            x[i] = (b[i] - sum) / diag;
        }
    }
}

// --- Динамический метод Якоби ---
static int jacobiSolveDyn(const double* A, const double* b, double* x, int n,
                          double eps, int maxIter) {
    double* xOld = new double[n];

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
                    sum += A[idx(i, j, n)] * xOld[j];
                }
            }
            x[i] = (b[i] - sum) / A[idx(i, i, n)];
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

    delete [] xOld;
    return iter;
}

// --- Динамический метод Гаусса–Зейделя ---
static int gaussSeidelSolveDyn(const double* A, const double* b, double* x, int n,
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
                    sum += A[idx(i, j, n)] * x[j];
                }
            }
            double newXi = (b[i] - sum) / A[idx(i, i, n)];
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

// --- Основная функция бенчмарка ---
void runComplexityBenchmark() {
    ensureCsvDirBenchmark();

    ofstream file("csv/complexity.csv");
    if (!file.is_open()) {
        cout << "Не удалось открыть csv/complexity.csv для записи.\n";
        return;
    }

    file << "n,method,time_ms,iterations\n";

    // Размеры систем, по которым будем оценивать сложность
    const int sizes[] = {10, 20, 40, 80, 120};
    const int numSizes = sizeof(sizes) / sizeof(sizes[0]);

    double eps = 1e-3;
    int maxIter = 1000;

    srand((unsigned int)time(nullptr));

    for (int s = 0; s < numSizes; ++s) {
        int n = sizes[s];

        // Базовая система
        double* A0 = new double[n * n];
        double* b0 = new double[n];
        generateDiagonalDominantSystem(A0, b0, n);

        // Рабочие массивы
        double* A = new double[n * n];
        double* b = new double[n];
        double* x = new double[n];

        // ----- Метод Гаусса -----
        for (int i = 0; i < n * n; ++i) A[i] = A0[i];
        for (int i = 0; i < n; ++i) b[i] = b0[i];

        auto startG = chrono::high_resolution_clock::now();
        gaussSolveDyn(A, b, x, n);
        auto endG = chrono::high_resolution_clock::now();
        double timeG = chrono::duration<double, milli>(endG - startG).count();

        file << n << ",Gauss," << timeG << "," << 0 << "\n";

        // ----- Метод Якоби -----
        for (int i = 0; i < n * n; ++i) A[i] = A0[i];
        for (int i = 0; i < n;     ++i) b[i] = b0[i];

        auto startJ = chrono::high_resolution_clock::now();
        int iterJ = jacobiSolveDyn(A, b, x, n, eps, maxIter);
        auto endJ = chrono::high_resolution_clock::now();
        double timeJ = chrono::duration<double, milli>(endJ - startJ).count();

        file << n << ",Jacobi," << timeJ << "," << iterJ << "\n";

        // ----- Метод Гаусса–Зейделя -----
        for (int i = 0; i < n * n; ++i) A[i] = A0[i];
        for (int i = 0; i < n;     ++i) b[i] = b0[i];

        auto startGS = chrono::high_resolution_clock::now();
        int iterGS = gaussSeidelSolveDyn(A, b, x, n, eps, maxIter);
        auto endGS = chrono::high_resolution_clock::now();
        double timeGS = chrono::duration<double, milli>(endGS - startGS).count();

        file << n << ",Gauss-Seidel," << timeGS << "," << iterGS << "\n";

        delete [] A0;
        delete [] b0;
        delete [] A;
        delete [] b;
        delete [] x;
    }

    file.close();
    cout << "Бенчмарк завершён, результаты в csv/complexity.csv\n";
}
