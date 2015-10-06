'''
Este script resuelve el oscilador de van der Pol
usando un metodo de Runge-Kutta orden 3. Las
funciones implementadas estan directamente relacionadas
con la solucion de la ecuacion diferencial: f es el
lado derecho de la EDO, y el resto son funciones
propias del metodo Runge-Kutta orden 3.

Hacia el final del script se resuelve el oscilador
de van der Pol para dos condiciones iniciales
diferentes, graficando los resultados para y(s) y
para el espacio y,dy/ds. Se guardan las imagenes
en formato .eps.
'''

import numpy as np
import matplotlib.pyplot as plt

def f(y, dy, mu):
    return dy, -y -mu*(y**2 - 1)*dy

def get_K1(y_n, dy_n, mu, h, f):
    f_eval=f(y_n, dy_n, mu)
    return h*f_eval[0], h*f_eval[1]

def get_K2(y_n, dy_n, mu, h, f, K1):
    f_eval=f(y_n + K1[0]/2., dy_n + K1[1]/2., mu)
    return h*f_eval[0], h*f_eval[1]

def get_K3(y_n, dy_n, mu, h, f, K1, K2):
    f_eval=f(y_n - K1[0] - 2*K2[0], dy_n + - K1[1] - 2*K2[1], mu)
    return h*f_eval[0], h*f_eval[1]

def RK3_step(y_n, dy_n, mu, h, f):
    K1=get_K1(y_n, dy_n, mu, h, f)
    K2=get_K2(y_n, dy_n, mu, h, f, K1)
    K3=get_K3(y_n, dy_n, mu, h, f, K1, K2)
    y_n1 = y_n + (K1[0] + 4*K2[0] + K3[0])/6.
    dy_n1 = dy_n + (K1[1] + 4*K2[1] + K3[1])/6.
    return y_n1, dy_n1

mu = 1.350
N_steps = 40000
h = 20*np.pi/N_steps
s = np.linspace(0,20*np.pi, N_steps+1)
y = np.zeros(N_steps+1)
dy = np.zeros(N_steps+1)

##Caso 1
y[0]=0.1
dy[0]=0.
for i in range(N_steps):
    y[i+1], dy[i+1] = RK3_step(y[i], dy[i], mu, h, f)

plt.figure(1)
plt.plot(s,y)
plt.xlabel('$s$')
plt.ylabel('$y(s)$')
plt.title('$y(s)$ para $\mu^* = 1.350$, $y(s=0)=0.1$, $dy/ds(s=0)=0$')
plt.savefig('ys1.eps')

plt.figure(2)
plt.plot(y,dy)
plt.xlabel('$y(s)$')
plt.ylabel('$dy/ds(s)$')
plt.title('$y(s)$ para $\mu^* = 1.350$, $y(s=0)=0.1$, $dy/ds(s=0)=0$')
plt.savefig('ydy1.eps')

#Caso 2
y[0]=4.0
dy[0]=0.
for i in range(N_steps):
    y[i+1], dy[i+1] = RK3_step(y[i], dy[i], mu, h, f)

plt.figure(3)
plt.plot(s,y)
plt.xlabel('$s$')
plt.ylabel('$y(s)$')
plt.title('$y(s)$ para $\mu^* = 1.350$, $y(s=0)=4.0$, $dy/ds(s=0)=0$')
plt.savefig('ys2.eps')

plt.figure(4)
plt.plot(y,dy)
plt.xlabel('$y(s)$')
plt.ylabel('$dy/ds(s)$')
plt.title('$y(s)$ para $\mu^* = 1.350$, $y(s=0)=4.0$, $dy/ds(s=0)=0$')
plt.savefig('ydy2.eps')

plt.show()
plt.draw()
