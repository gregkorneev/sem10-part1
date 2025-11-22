import csv
import os


def read_all_methods(path: str):
    # Определяем разделитель по первой строке
    with open(path, "r") as f:
        first_line = f.readline()
        delimiter = ";" if ";" in first_line and "," not in first_line else ","

    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            if not row:
                continue

            method = (row.get("method") or row.get("Method") or "").strip()
            if not method:
                continue

            it_str = row.get("iterations") or row.get("Iterations") or "0"
            res_str = row.get("residual") or row.get("Residual") or "0"
            err_str = row.get("error_vs_gauss") or row.get("Error_vs_gauss") or "0"

            try:
                iterations = int(it_str)
                residual = float(res_str.replace(",", "."))
                error_vs_gauss = float(err_str.replace(",", "."))
            except ValueError:
                # если вдруг попалась “грязная” строка — пропускаем
                continue

            rows.append((method, iterations, residual, error_vs_gauss))
    return rows


def main():
    csv_path = os.path.join("csv", "all_methods.csv")

    if not os.path.exists(csv_path):
        print(f"Файл {csv_path} не найден")
        return

    rows = read_all_methods(csv_path)
    if not rows:
        print("Не удалось прочитать данные из csv/all_methods.csv")
        return

    print("Оценка методов (по данным all_methods.csv):\n")
    print(f"{'Метод':15s} {'Итерации':>10s} {'Невязка':>15s} {'Погрешность':>15s}")
    print("-" * 60)
    for method, iters, residual, error in rows:
        print(f"{method:15s} {iters:10d} {residual:15.3e} {error:15.3e}")


if __name__ == "__main__":
    main()
