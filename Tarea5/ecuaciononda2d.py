import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

Nt = 400
Nx = 20
Ny = 20

x_vals = np.linspace(0, 2, Nx)
y_vals = np.linspace(0, 2, Ny)
time_vals = np.linspace(0, 3, Nt)

delta_x = x_vals[1] - x_vals[0]
delta_y = y_vals[1] - y_vals[0]
delta_t = time_vals[1] - time_vals[0]

velocity = 2.

lambda_ = velocity * delta_t / delta_x
mu_ = velocity * delta_t / delta_y

print(lambda_, mu_)

def initial_condition(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)

solution = np.zeros((Nt, Nx, Ny))
for i in range(len(x_vals)):
    for j in range(len(y_vals)):
        solution[0, i, j] = initial_condition(x_vals[i], y_vals[j])

def compute_solution():
    damping_factor = 0 * delta_t

    for l in range(1, len(time_vals)):
        if l == 1:
            solution[l, :, :] = solution[l - 1, :, :]
        else:
            for i in range(1, len(x_vals) - 1):
                for j in range(1, len(y_vals) - 1):
                    solution[l, i, j] = 2 * (1 - lambda_ ** 2 - mu_ ** 2) * solution[l - 1, i, j] \
                                + lambda_ ** 2 * (solution[l - 1, i + 1, j] + solution[l - 1, i - 1, j]) \
                                + mu_ ** 2 * (solution[l - 1, i, j + 1] + solution[l - 1, i, j - 1]) \
                                - solution[l - 2, i, j] \
                                - damping_factor * solution[l - 1, i, j] + damping_factor * solution[l - 2, i, j]

compute_solution()

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(x_vals, y_vals)

def initialize_plot():
    ax.set_xlim3d(0, 2)
    ax.set_ylim3d(0, 2)
    ax.set_zlim3d(-1, 1)

def update_plot(i):
    ax.clear()
    initialize_plot()
    ax.plot_surface(X, Y, solution[i, :, :], cmap='viridis')

animation_obj = animation.FuncAnimation(fig, update_plot, frames=len(time_vals), init_func=initialize_plot)
plt.show()
