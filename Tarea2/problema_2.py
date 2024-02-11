import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(u, t, q):
    return u**q

q_values = [0., 0.2, 0.4, 0.7, 0.9, 1.]


u0 = 1.0
t = np.linspace(0, 10, 100)  # Rango de tiempo de 0 a 10
for q in q_values:
    u = odeint(model, u0, t, args=(q,))
    plt.plot(t, u, label=f'q = {q}')

# Configurar el gráfico
plt.title('Solución numérica de du/dt = u^q')
plt.xlabel('Tiempo')
plt.ylabel('u(t)')
plt.legend()
plt.grid(True)
plt.show()
