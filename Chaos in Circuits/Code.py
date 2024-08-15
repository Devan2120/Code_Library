import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# parameters
alpha = 15.395
beta = 28
R = -1.143
C_2 = -0.714

def chua(u, t):
    x, y, z = u  
    # electrical response of the nonlinear resistor
    f_x = C_2*x + 0.5*(R-C_2)*(abs(x+1)-abs(x-1))
    dudt = [alpha*(y-x-f_x), x - y + z, -beta * y]
    return dudt

# time discretization
t_0 = 0
dt = 1e-3
t_final = 300
t = np.arange(t_0, t_final, dt)

# initial conditions
u0 = [0.1,0,0]

# integrate ode system
sol = odeint(chua, u0, t)

# 3d-plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.plot(sol[:,0],sol[:,1],sol[:,2])
plt.show()