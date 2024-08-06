import numpy as np
import matplotlib.pyplot as plt

def read_data(filename):
    data = np.loadtxt(filename, skiprows=1)  
    t = data[:, 0]
    y_rk4 = data[:, 1]
    y_exact = data[:, 2]
    return t, y_rk4, y_exact

def calculate_error(y_rk4, y_exact):
    return np.abs(y_rk4 - y_exact)

def plot_results(t, y_rk4, y_exact, error):
    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.plot(t, y_rk4, 'b--', label='RK4')
    plt.plot(t, y_exact, 'r-', label='Exact')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.legend()
    plt.title('RK4 vs Exact Solution')

    plt.subplot(2, 1, 2)
    plt.plot(t, error, 'k-')
    plt.xlabel('t')
    plt.ylabel('Error')
    plt.title('Error between RK4 and Exact Solution')

    plt.tight_layout()
    plt.show()

def main():
    filename = 'output_h_0.100.txt'
    t, y_rk4, y_exact = read_data(filename)
    error = calculate_error(y_rk4, y_exact)
    plot_results(t, y_rk4, y_exact, error)

if __name__ == '__main__':
    main()
