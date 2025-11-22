import os
import matplotlib.pyplot as plt

def main():
    # Размеры матрицы
    n_values = [10, 20, 50, 100, 200, 500]

    # Число операций (условные единицы)
    # Гаусс ~ O(n^3), итерационные ~ O(n^2)
    gauss_ops = [0.5 * n**3 for n in n_values]
    iter_ops  = [5.0 * n**2 for n in n_values]

    os.makedirs("plots", exist_ok=True)

    plt.figure(figsize=(10, 4))

    plt.semilogy(n_values, gauss_ops, marker="o", label="Метод Гаусса O(n³)")
    plt.semilogy(n_values, iter_ops,  marker="o", label="Итерационные методы O(n²)")

    plt.title("Сравнение вычислительной сложности методов")
    plt.xlabel("Размер матрицы n")
    plt.ylabel("Количество операций (условные единицы)")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()

    plt.tight_layout()
    plt.savefig("plots/complexity_theory.png", dpi=200)
    plt.close()

if __name__ == "__main__":
    main()
