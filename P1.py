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
plt.figure(2)
plt.plot(y,dy)


#Caso 2
y[0]=4.0
dy[0]=0.
for i in range(N_steps):
    y[i+1], dy[i+1] = RK3_step(y[i], dy[i], mu, h, f)

plt.figure(3)
plt.plot(s,y)
plt.figure(4)
plt.plot(y,dy)

plt.show()
plt.draw()
