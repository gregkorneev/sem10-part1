#!/usr/bin/env python3
import csv
import os
import math

import matplotlib.pyplot as plt


def read_complexity_csv(path: str):
    """
    Читаем csv/complexity.csv формата:
    n,method,time_ms,iterations
    """
    data = {}  # method -> {"n": [...], "time": [...]}

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            method = row["method"]
            n = int(row["n"])
            t = float(row["time_ms"])

            if method not in data:
                data[method] = {"n": [], "time": []}
            data[method]["n"].append(n)
            data[method]["time"].append(t)

    return data


def ensure_plots_dir():
    os.makedirs("plots", exist_ok=True)


def estimate_power_trend(ns, times):
    """
    Оцениваем степенной тренд T(n) ≈ C * n^k.
    Делается линейная регрессия в лог-лог координатах:
    ln T = k ln n + ln C.
    Возвращаем (C, k).
    """
    xs = [math.log(n) for n in ns]
    ys = [math.log(t) for t in times]

    m = len(xs)
    if m < 2:
        return None, None

    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))

    denom = m * sxx - sx * sx
    if abs(denom) < 1e-12:
        return None, None

    k = (m * sxy - sx * sy) / denom       # степень
    b = (sy - k * sx) / m                 # ln C
    C = math.exp(b)

    return C, k


def plot_all_together(data):
    """
    Общий график (как раньше, без трендов).
    """
    plt.figure()

    for method, vals in data.items():
        ns = vals["n"]
        times = vals["time"]
        plt.plot(ns, times, marker="o", label=method)

    plt.title("Зависимость времени работы от размера системы")
    plt.xlabel("Размер системы n")
    plt.ylabel("Время, мс")

    plt.xscale("log")
    plt.yscale("log")

    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/complexity_time.png")
    plt.close()


def plot_power_trends(data, trends):
    """
    Для каждого метода: отдельный график с точками и степенным трендом
    T(n) = C * n^k, формула выводится на графике.
    """
    for method, vals in data.items():
        ns = vals["n"]
        times = vals["time"]
        C, k = trends[method]

        plt.figure()

        # точки бенчмарка
        plt.scatter(ns, times, label="измерения")

        # степенной тренд
        if C is not None and k is not None:
            # строим гладкую линию по тем же n
            ns_sorted = sorted(ns)
            fit_times = [C * (n ** k) for n in ns_sorted]
            plt.plot(ns_sorted, fit_times, label="степенной тренд")

        plt.title(f"Время работы метода {method}")
        plt.xlabel("Размер системы n")
        plt.ylabel("Время, мс")

        plt.xscale("log")
        plt.yscale("log")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)

        # подпись формулы
        if C is not None and k is not None:
            # красиво округлим C, чтобы формула была читаемой
            if C >= 1:
                C_str = f"{C:.3g}"
            else:
                C_str = f"{C:.2e}"
            text = f"T(n) ≈ {C_str} · n^{k:.2f}\n(время в мс)"
            ax = plt.gca()
            ax.text(
                0.05, 0.95,
                text,
                transform=ax.transAxes,
                va="top",
                ha="left",
            )

        plt.legend()
        plt.tight_layout()

        safe_name = method.replace(" ", "_")
        out_path = os.path.join("plots", f"complexity_power_{safe_name}.png")
        plt.savefig(out_path)
        plt.close()


def write_md_description(trends):
    """
    MD-файл с формулами степенного тренда.
    """
    md_path = os.path.join("plots", "complexity_time.md")

    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# Степенные тренды сложности методов\n\n")
        f.write(
            "Для каждого метода по результатам бенчмарка оценена зависимость\n"
            "вида `T(n) ≈ C · n^k`, где `T` — время в миллисекундах, `n` — размер системы.\n\n"
        )

        for method, (C, k) in trends.items():
            if C is None or k is None:
                continue
            if C >= 1:
                C_str = f"{C:.3g}"
            else:
                C_str = f"{C:.2e}"
            f.write(f"- **{method}**: T(n) ≈ {C_str} · n^{k:.2f}\n")

        f.write(
            "\nКоэффициент k показывает степень роста времени от n.\n"
            "Для метода Гаусса k близок к 3 → сложность порядка O(n³).\n"
            "Для методов Якоби и Гаусса–Зейделя k близок к 2 → сложность порядка O(k · n²).\n"
        )


def main():
    csv_path = os.path.join("csv", "complexity.csv")
    if not os.path.exists(csv_path):
        print("Файл csv/complexity.csv не найден.")
        print("Сначала запусти ./build/app, чтобы программа создала CSV.")
        return

    data = read_complexity_csv(csv_path)
    ensure_plots_dir()

    # оценка степенного тренда для каждого метода
    trends = {}
    print("Степенные тренды T(n) ≈ C · n^k (T в мс):")
    for method, vals in data.items():
        C, k = estimate_power_trend(vals["n"], vals["time"])
        trends[method] = (C, k)
        if C is not None and k is not None:
            if C >= 1:
                C_str = f"{C:.3g}"
            else:
                C_str = f"{C:.2e}"
            print(f"{method}: T(n) ≈ {C_str} · n^{k:.2f}")
        else:
            print(f"{method}: недостаточно данных для оценки.")

    # общий график
    plot_all_together(data)

    # отдельные графики с степенным трендом
    plot_power_trends(data, trends)

    # md-описание
    write_md_description(trends)

    print("Графики сохранены в папке plots/:")
    print("  complexity_time.png (общий)")
    for method in data.keys():
        safe_name = method.replace(' ', '_')
        print(f"  complexity_power_{safe_name}.png (степенной тренд)")
    print("Формулы тренда записаны в plots/complexity_time.md")


if __name__ == "__main__":
    main()
