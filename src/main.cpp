#include <iostream>
#include <cmath>
#include <iomanip>
#include "../include/methods.h"
#include "../include/csv_output.h"
#include "../include/benchmark.h"

using namespace std;

int main() {
    setlocale(LC_ALL, "Russian");

    // Система 14:
    // 4x - y + 2z = 11
    // x + 3y - z = 4
    // 2x + y + z = 7
    double A[N][N] = {
        {4.0, -1.0,  2.0},
        {1.0,  3.0, -1.0},
        {2.0,  1.0,  1.0}
    };

    // Вектор правых частей для системы 14
    double b[N] = {11.0, 4.0, 7.0};

    // ---------- 1. Метод Гаусса ----------
    double xGauss[N];

    cout << "=== Метод Гаусса ===\n";

    if (!gaussSolve(A, b, xGauss, N)) {
        cout << "Решить систему методом Гаусса не удалось.\n";
        return 1;
    }

    cout << "Решение (Гаусс):\n";
    printVector(xGauss, N);

    double resGauss = residualInfinityNorm(A, xGauss, b, N);
    cout << "Невязка (норма max) для метода Гаусса: " << resGauss << "\n\n";

    // ---------- 2. Метод Якоби ----------
    cout << "=== Метод Якоби ===\n";

    double xJacobi[N];
    double eps = 1e-3;
    int maxIter = 1000; // чтобы точно успеть сойтись

    int iterJacobi = jacobiSolve(A, b, xJacobi, N, eps, maxIter);

    cout << "Число итераций: " << iterJacobi << "\n";
    cout << "Решение (Якоби):\n";
    printVector(xJacobi, N);

    double resJacobi = residualInfinityNorm(A, xJacobi, b, N);
    cout << "Невязка (норма max): " << resJacobi << "\n";

    double errJacobi = errorNorm2(xJacobi, xGauss, N);
    cout << "Погрешность относительно решения метода Гаусса (норма 2): "
         << errJacobi << "\n\n";

    // ---------- 3. Метод Гаусса–Зейделя ----------
    cout << "=== Метод Гаусса–Зейделя ===\n";

    double xGS[N];
    int iterGS = gaussSeidelSolve(A, b, xGS, N, eps, maxIter);

    cout << "Число итераций: " << iterGS << "\n";
    cout << "Решение (Гаусса–Зейделя):\n";
    printVector(xGS, N);

    double resGS = residualInfinityNorm(A, xGS, b, N);
    cout << "Невязка (норма max): " << resGS << "\n";

    double errGS = errorNorm2(xGS, xGauss, N);
    cout << "Погрешность относительно решения метода Гаусса (норма 2): "
         << errGS << "\n\n";

    cout << "=== Оценка вычислительной сложности ===\n";
    cout << "Метод Гаусса: O(n^3) для n x n.\n";
    cout << "Метод Якоби: ~ O(k * n^2), где k - число итераций.\n";
    cout << "Метод Гаусса–Зейделя: ~ O(k * n^2).\n\n";

    // --- Запись результатов в CSV ---
    writeAllMethodsCsv(
        xGauss, resGauss,
        xJacobi, iterJacobi, resJacobi, errJacobi,
        xGS, iterGS, resGS, errGS
    );

    writeGaussCsv(xGauss, resGauss);
    writeJacobiCsv(xJacobi, iterJacobi, resJacobi, errJacobi);
    writeGaussSeidelCsv(xGS, iterGS, resGS, errGS);

    // --- Бенчмарк сложности (время vs размер n) ---
    runComplexityBenchmark();

    return 0;
}
