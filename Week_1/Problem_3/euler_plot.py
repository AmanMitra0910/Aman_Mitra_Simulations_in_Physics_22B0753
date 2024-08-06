import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    data = np.loadtxt(filename , skiprows = 1)
    t = data[:, 0]
    y_euler = data[:, 1]
    y_exact = data[:, 2]
    return t , y_euler , y_exact

t_01, y_euler_01, y_exact_01 = read_data('output_h_0.100.txt')
t_001, y_euler_001, y_exact_001 = read_data('output_h_0.010.txt')

plt.figure(figsize = (12,6))

plt.plot(t_01, y_euler_01, 'o-', label='Euler h=0.1', markersize=4)
plt.plot(t_001, y_euler_001, 's-', label='Euler h=0.o1', markersize=2)
plt.plot(t_001, y_exact_001, '-', label='Exact Solution', linewidth=2)

plt.xlabel('Time t')
plt.ylabel('y(t)')
plt.title('Comparison of Euler Method and Exact Solution')
plt.legend()
plt.grid()
plt.show()