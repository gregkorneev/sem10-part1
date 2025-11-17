#ifndef CSV_OUTPUT_H
#define CSV_OUTPUT_H

void writeAllMethodsCsv(
    const double xGauss[], double resGauss,
    const double xJacobi[], int iterJacobi, double resJacobi, double errJacobi,
    const double xGS[], int iterGS, double resGS, double errGS
);

void writeGaussCsv(const double x[], double residual);
void writeJacobiCsv(const double x[], int iterations, double residual, double error);
void writeGaussSeidelCsv(const double x[], int iterations, double residual, double error);

#endif // CSV_OUTPUT_H
