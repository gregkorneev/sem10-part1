# Семинар 10 — Решение СЛАУ тремя методами (Гаусс, Якоби, Гаусса–Зейделя)

```bash
rm -rf build && rm -rf csv && rm -rf plots  && cmake -B build && cmake --build build && ./build/app && python3 .py/plot_results.py && python3 .py/plot_complexity.py && python3 .py/plot_complexity_time_log.py && python3 .py/plot_complexity_theory.py && python3 .py/print_method_stats.py
