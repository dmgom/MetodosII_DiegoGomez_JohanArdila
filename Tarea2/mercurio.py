import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Planeta:
    def __init__(self, e, a, dt, alpha):
        self.dt = dt
        self.e = e  # Excentricidad
        self.a = a  # Semi-eje mayor
        self.alpha = alpha  # Valor de alpha

        self.G = 4 * np.pi**2  # Unidades gaussianas

        self.r = np.zeros(2)
        self.v = np.zeros_like(self.r)

        # Inicializar en el afelio
        self.r[0] = self.a * (1 + self.e)
        self.v[1] = np.sqrt(self.G * (1 - self.e) / (self.a * (1 + self.e)))

        # Registro de tiempo y ángulo
        self.time_to_perihelion = []
        self.arrival_angle = []

    def acceleration(self, r):
        r_norm = np.linalg.norm(r)
        return -self.G * r / r_norm**3 * (1 + self.alpha / r_norm**2)

    def evolve(self):
        # Calcular la aceleración en el paso actual
        a_current = self.acceleration(self.r)
        
        # Actualizar la posición usando el método de Verlet modificado
        self.r += self.v * self.dt + 0.5 * a_current * self.dt**2
        
        # Calcular la aceleración en el paso siguiente
        a_next = self.acceleration(self.r)
        
        # Actualizar la velocidad usando el método de Verlet modificado
        self.v += 0.5 * (a_current + a_next) * self.dt

        # Verificar si llegamos al perihelio
        if np.isclose(self.r[0], self.a * (1 - self.e), rtol=1e-4):
            # Calcular el ángulo de llegada
            angle = np.arctan2(self.r[1], self.r[0]) * 180 / np.pi
            self.arrival_angle.append(angle)
            print(arrival_angle)
            # Registrar el tiempo
            self.time_to_perihelion.append(len(self.time_to_perihelion) * self.dt)

def simulate_mercury(num_steps, dt, alpha):
    # Parámetros orbitales de Mercurio
    e_mercury = 0.205630
    a_mercury = 0.387098
    

    mercury = Planeta(e_mercury, a_mercury, dt, alpha)
   
    positions = []
    
    for _ in range(num_steps):
        mercury.evolve()
        positions.append(np.copy(mercury.r))
    return np.array(positions), mercury.time_to_perihelion, mercury.arrival_angle

num_steps = 10000
dt = 1e-4  # Paso temporal del mismo orden de alpha
alpha = 1.1e-8  # Valor de alpha


positions, time_to_perihelion, arrival_angle = simulate_mercury(num_steps, dt, alpha)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))


ax1.set_xlim(-0.5, 0.5)
ax1.set_ylim(-0.5, 0.5)
ax1.set_aspect('equal', adjustable='box')

line, = ax1.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,


def animate(i):
    line.set_data(positions[:i, 0], positions[:i, 1])
    return line,


ani = animation.FuncAnimation(fig, animate, init_func=init, frames=num_steps, interval=10, blit=True)

plt.show()
