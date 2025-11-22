#!/usr/bin/env python3
import csv
import os

import matplotlib.pyplot as plt


def read_all_methods_csv(path: str):
    methods = []
    iterations = []
    residuals = []
    errors = []

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            methods.append(row["method"])
            iterations.append(int(row["iterations"]))
            residuals.append(float(row["residual"]))
            errors.append(float(row["error_vs_gauss"]))

    return methods, iterations, residuals, errors


def ensure_plots_dir():
    os.makedirs("plots", exist_ok=True)


def plot_iterations(methods, iterations):
    plt.figure()
    plt.bar(methods, iterations)
    plt.title("Количество итераций по методам")
    plt.xlabel("Метод")
    plt.ylabel("Число итераций")
    plt.tight_layout()
    plt.savefig("plots/iterations.png")
    plt.close()


def plot_residuals(methods, residuals):
    plt.figure()
    plt.bar(methods, residuals)
    plt.title("Невязка (норма max) по методам")
    plt.xlabel("Метод")
    plt.ylabel("Невязка")
    plt.yscale("log")  # логарифмическая шкала, чтобы различия было видно
    plt.tight_layout()
    plt.savefig("plots/residuals.png")
    plt.close()


def plot_errors(methods, errors):
    plt.figure()
    plt.bar(methods, errors)
    plt.title("Погрешность относительно метода Гаусса")
    plt.xlabel("Метод")
    plt.ylabel("Погрешность (норма 2)")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig("plots/errors.png")
    plt.close()


def main():
    csv_path = os.path.join("csv", "all_methods.csv")

    if not os.path.exists(csv_path):
        print(f"Файл {csv_path} не найден. Сначала запусти ./build/app, "
              f"чтобы сгенерировать CSV.")
        return

    methods, iterations, residuals, errors = read_all_methods_csv(csv_path)
    ensure_plots_dir()

    plot_iterations(methods, iterations)
    plot_residuals(methods, residuals)
    plot_errors(methods, errors)

    print("Графики сохранены в папку plots/:")
    print("  plots/iterations.png")
    print("  plots/residuals.png")
    print("  plots/errors.png")


if __name__ == "__main__":
    main()
