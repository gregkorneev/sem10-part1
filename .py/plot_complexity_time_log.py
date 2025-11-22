import csv
import os
import math
import matplotlib.pyplot as plt


def read_complexity_csv(path):
    """
    Ожидаемый формат csv/complexity.csv:

    N,Method,Time_ms,Iterations
    10,Gauss,0.01,0
    10,Jacobi,0.01,5
    10,Gauss-Seidel,0.0,4
    ...

    N        – размер матрицы
    Method   – 'Gauss', 'Jacobi', 'Gauss-Seidel'
    Time_ms  – время выполнения в миллисекундах
    """

    # Определяем разделитель
    with open(path, "r") as f:
        first_line = f.readline()
        delimiter = ";" if ";" in first_line and "," not in first_line else ","

    # словари: method -> список (N, time)
    data = {
        "Gauss": [],
        "Jacobi": [],
        "Gauss-Seidel": [],
    }

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            if not row:
                continue

            # названия столбцов точно такие, как в файле
            n_str = row.get("N") or row.get("n")
            method = (row.get("Method") or row.get("method") or "").strip()
            time_str = row.get("Time_ms") or row.get("time_ms") or row.get("Time") or row.get("time")

            if not (n_str and method and time_str):
                continue

            if method not in data:
                # если вдруг там какая-то другая строка – пропустим
                continue

            try:
                n = int(n_str)
                t = float(time_str.replace(",", "."))
            except ValueError:
                continue

            # для логарифма нельзя время <= 0
            if t <= 0.0:
                t = 1e-6

            data[method].append((n, t))

    # сортируем по N и раскладываем в списки
    def prepare_series(method_name):
        series = sorted(data[method_name], key=lambda x: x[0])
        n_vals = [n for n, _ in series]
        times = [t for _, t in series]
        return n_vals, times

    n_gauss, t_gauss = prepare_series("Gauss")
    n_jacobi, t_jacobi = prepare_series("Jacobi")
    n_gs, t_gs = prepare_series("Gauss-Seidel")

    # считаем, что N одни и те же, берём из Gauss
    return n_gauss, t_gauss, t_jacobi, t_gs


def plot_log_log_time(n_vals, gauss, jacobi, gs, out_path):
    """Строит график в логарифмических координатах по времени выполнения."""
    log_n = [math.log10(x) for x in n_vals]
    log_gauss = [math.log10(t) for t in gauss]
    log_jacobi = [math.log10(t) for t in jacobi]
    log_gs = [math.log10(t) for t in gs]

    plt.figure(figsize=(10, 4))

    plt.plot(log_n, log_gauss, marker="o", label="Метод Гаусса")
    plt.plot(log_n, log_jacobi, marker="o", label="Метод Якоби")
    plt.plot(log_n, log_gs, marker="o", label="Метод Гаусса-Зейделя")

    plt.title("Логарифмическая шкала времени выполнения")
    plt.xlabel("log₁₀(размера матрицы n)")
    plt.ylabel("log₁₀(времени выполнения, мс)")

    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    csv_path = os.path.join("csv", "complexity.csv")
    plots_dir = "plots"
    os.makedirs(plots_dir, exist_ok=True)

    n_vals, t_gauss, t_jacobi, t_gs = read_complexity_csv(csv_path)

    if not n_vals:
        print("⚠ Не удалось прочитать данные из csv/complexity.csv")
        return

    out_path = os.path.join(plots_dir, "complexity_time_log.png")
    plot_log_log_time(n_vals, t_gauss, t_jacobi, t_gs, out_path)
    print("✅ График сохранён в", out_path)


if __name__ == "__main__":
    main()
