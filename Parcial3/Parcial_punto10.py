import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
ν = 0.3
T = 10
Nx = Ny = 60
Nt = 500
L = 10
dx = dy = L / (Nx - 1)
dt = T / Nt
u = np.zeros((Nt, Nx, Ny))

x = np.linspace(-5, 5, Nx)
y = np.linspace(-5, 5, Ny)
X, Y = np.meshgrid(x, y)
u[0] = 5 * np.exp(-(X**2 + Y**2))
for n in range(Nt - 1):
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            u[n + 1, i, j] = (u[n, i, j] - 
                              u[n, i, j] * (dt / dx) * (u[n, i, j] - u[n, i - 1, j]) -
                              u[n, i, j] * (dt / dy) * (u[n, i, j] - u[n, i, j - 1]) +
                              ν * (dt / dx**2) * (u[n, i + 1, j] - 2 * u[n, i, j] + u[n, i - 1, j]) +
                              ν * (dt / dy**2) * (u[n, i, j + 1] - 2 * u[n, i, j] + u[n, i, j - 1]))
def animate(n):
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, u[n], cmap='viridis')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')
    ax.set_xlim3d(-4,4)
    ax.set_ylim3d(-4,4)
    ax.set_zlim3d(0,2)
fig = plt.figure()
ani = FuncAnimation(fig, animate, frames=Nt)
plt.show()
