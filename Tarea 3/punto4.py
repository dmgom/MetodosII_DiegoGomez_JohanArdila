import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
G = 6.67e-11 
m_T = 5.9736e24 
r_T = 6.3781e6 
m_L = 0.07349e24 
r_L = 1.7374e6 
d = 3.844e8 
omega = 2.6617e-6 
Delta = G * m_T / d**3 
mu = m_L / m_T 

def f(t, y):
    r, phi, pr, pphi = y
    r_prime = np.sqrt(1 + r**2 - 2*r*np.cos(phi - omega*t))
    drdt = pr
    dphidt = pphi / r**2
    dpdr = pphi**2 / r**3 - Delta / r**2 + mu / r_prime**3 * (r - np.cos(phi - omega*t))
    dpdphi = -Delta * mu * r / r_prime**3 * np.sin(phi - omega*t)
    return [drdt, dphidt, dpdr, dpdphi]

def RK4(f, t0, tf, y0, h):
    t = t0
    y = y0
    sol = [[t, *y]]
    while t < tf:
        k1 = np.array(f(t, y)) * h
        k2 = np.array(f(t + 0.5*h, y + 0.5*k1)) * h
        k3 = np.array(f(t + 0.5*h, y + 0.5*k2)) * h
        k4 = np.array(f(t + h, y + k3)) * h
        y += (k1 + 2*k2 + 2*k3 + k4) / 6
        t += h
        sol.append([t, *y])
    return np.array(sol)

v_0 = np.sqrt(G * m_T / d)  # Velocidad inicial cercana a la velocidad de escape de la Tierra
theta = np.pi / 4  # Ángulo de lanzamiento (rad)
phi_0 = np.pi / 4  # Latitud inicial (rad)
p_r0 = v_0 * np.cos(theta - phi_0)
p_phi0 = (r_T / d)**2
y0 = [r_T / d, phi_0, p_r0, p_phi0]
t0 = 0
tf = 10*24*3600  # 10 días en segundos
h = 1000  # Paso de tiempo (segundos)

sol = RK4(f, t0, tf, y0, h)
def init():
    line.set_data([], [])
    return line,

# Función de animación
def animate(i):
    x = sol[:i, 1] * d * np.cos(sol[:i, 2])  
    y = sol[:i, 1] * d * np.sin(sol[:i, 2])  
    line.set_data(x, y)
    return line,

# Configurar la figura y la animación
fig, ax = plt.subplots()
ax.set_xlim(-d, d)
ax.set_ylim(-d, d)
ax.set_xlabel('Posición x (m)')
ax.set_ylabel('Posición y (m)')
ax.set_title('Trayectoria de la nave espacial hacia la Luna')
line, = ax.plot([], [], lw=2)
ani = FuncAnimation(fig, animate, init_func=init, frames=len(sol), blit=True)

# Mostrar la animación
plt.show()