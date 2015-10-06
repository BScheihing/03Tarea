'''
Este script resuelve el sistema de Lorenz usando
el integrador dopri5 de la libreria scipy.integrate.
Los parametros asociados al sistema de Lorenz estan
fijos en la implementacion, para obtener la solucion
conocida como el Atractor de Lorenz.

La funcion f define el sistema con los parametros ya
fijados. Luego el script integra (x,y,z) desde t=0
hasta t=100 y crea una figura con la trayectoria en
3D. Al final del script se guarda una de sus vistas
en una imagen .eps.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

def f(t, xyz):
    dx = 10.0 * (xyz[1] - xyz[0])
    dy = xyz[0] * (28.0 - xyz[2]) - xyz[1]
    dz = xyz[0] * xyz[1] - 8.0/3.0 * xyz[2]
    return [dx, dy, dz]

tf = 100.
N_steps = 10000
dt = tf/N_steps
x=np.zeros(N_steps+1)
y=np.zeros(N_steps+1)
z=np.zeros(N_steps+1)
xyz0 = [1, 1, 1]
t=np.linspace(0, 10, N_steps+1)

n=1
solver = ode(f)
solver.set_integrator('dopri5', atol=1E-6, rtol=1E-4)
solver.set_initial_value(xyz0)
while solver.successful() and solver.t < tf and n<=N_steps:
    solver.integrate(solver.t+dt)
    t[n] = solver.t
    x[n] = solver.y[0]
    y[n] = solver.y[1]
    z[n] = solver.y[2]
    n+=1

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')
ax.plot(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Atractor de Lorenz, $(x_0,y_0,z_0) = (1,1,1)$')
plt.savefig('Lorenz.eps')

plt.show()
