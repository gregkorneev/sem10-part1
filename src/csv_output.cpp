#include <fstream>
#include <iostream>
#include <cerrno>

#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #include <sys/types.h>
    #define MKDIR(path) mkdir(path, 0755)
#endif

#include "../include/methods.h"


using namespace std;

// Создание папки csv (если её ещё нет)
void ensureCsvDir() {
    int result = MKDIR("csv");
    // Если директория уже существует, это не ошибка.
    // На Unix это errno == EEXIST, на Windows _mkdir тоже вернёт ошибку,
    // но для нашей простой лабы можно её тихо игнорировать.
    (void)result;
}

// Общая таблица по всем методам
void writeAllMethodsCsv(
    const double xGauss[], double resGauss,
    const double xJacobi[], int iterJacobi, double resJacobi, double errJacobi,
    const double xGS[], int iterGS, double resGS, double errGS
) {
    ensureCsvDir();

    ofstream file("csv/all_methods.csv");
    if (!file.is_open()) {
        cout << "Не удалось открыть файл csv/all_methods.csv для записи.\n";
        return;
    }

    file << "method,iterations,x1,x2,x3,residual,error_vs_gauss\n";

    file << "Gauss,"
         << 0 << ","
         << xGauss[0] << ","
         << xGauss[1] << ","
         << xGauss[2] << ","
         << resGauss << ","
         << 0.0 << "\n";

    file << "Jacobi,"
         << iterJacobi << ","
         << xJacobi[0] << ","
         << xJacobi[1] << ","
         << xJacobi[2] << ","
         << resJacobi << ","
         << errJacobi << "\n";

    file << "Gauss-Seidel,"
         << iterGS << ","
         << xGS[0] << ","
         << xGS[1] << ","
         << xGS[2] << ","
         << resGS << ","
         << errGS << "\n";

    file.close();
}

// Отдельный CSV для метода Гаусса
void writeGaussCsv(const double x[], double residual) {
    ensureCsvDir();

    ofstream file("csv/gauss.csv");
    if (!file.is_open()) {
        cout << "Не удалось открыть файл csv/gauss.csv для записи.\n";
        return;
    }

    file << "method,iterations,x1,x2,x3,residual,error_vs_gauss\n";
    file << "Gauss,"
         << 0 << ","
         << x[0] << ","
         << x[1] << ","
         << x[2] << ","
         << residual << ","
         << 0.0 << "\n";

    file.close();
}

// Отдельный CSV для метода Якоби
void writeJacobiCsv(const double x[], int iterations, double residual, double error) {
    ensureCsvDir();

    ofstream file("csv/jacobi.csv");
    if (!file.is_open()) {
        cout << "Не удалось открыть файл csv/jacobi.csv для записи.\n";
        return;
    }

    file << "method,iterations,x1,x2,x3,residual,error_vs_gauss\n";
    file << "Jacobi,"
         << iterations << ","
         << x[0] << ","
         << x[1] << ","
         << x[2] << ","
         << residual << ","
         << error << "\n";

    file.close();
}

// Отдельный CSV для метода Гаусса–Зейделя
void writeGaussSeidelCsv(const double x[], int iterations, double residual, double error) {
    ensureCsvDir();

    ofstream file("csv/gauss_seidel.csv");
    if (!file.is_open()) {
        cout << "Не удалось открыть файл csv/gauss_seidel.csv для записи.\n";
        return;
    }

    file << "method,iterations,x1,x2,x3,residual,error_vs_gauss\n";
    file << "Gauss-Seidel,"
         << iterations << ","
         << x[0] << ","
         << x[1] << ","
         << x[2] << ","
         << residual << ","
         << error << "\n";

    file.close();
}
