import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

n, dt = 100, 0.01
'''
R, alpha, beta = 4, 30, 60
x_dir, y_dir, z_dir = 1, 1, 1
u, l, m, w, zr, w1 = 0.5, 5, 1.45, 1, 1, 1
u0, u1, kt = 100, 0, 69



k = 0.6 * u0
stx, stp = (kt / 2 / k) ** 0.5, (m * kt) ** 0.5
n_x = np.random.normal(0, stx, n)
n_vx = np.random.normal(0, stp, n) / m
n_y = np.random.normal(0, stx, n)
n_vy = np.random.normal(0, stp, n) / m
n_z = np.random.normal(0, stx, n)
n_vz = np.random.normal(0, stp, n) / m
tau, count = 2 * l / u, 0
t = np.arange(0, tau, dt)

var = np.sqrt(x_dir **2 + y_dir **2 + z_dir ** 2)
cosa, cosb, cosg = x_dir / var, y_dir / var, z_dir / var
x1, y1, z1 = R * np.cos(beta / 180 * 3.14) * np.cos(alpha / 180 * 3.14), R * np.cos(beta / 180 * 3.14) * np.sin(alpha / 180 * 3.14), R * np.sin(beta / 180 * 3.14)
'''


def dsdt(t, S):
  x, vx, y, vy, z, vz = S
  f =  np.exp(-(x **2 + y ** 2) / w / w - z ** 2 / zr / zr)
  deltax, deltay ,deltaz = -x1 - cosa * (u * t - l) + x, -y1 - cosb * (u * t - l) + y, -z1 - cosg * (u * t - l) + z
  exp = np.exp(-(np.sqrt(deltax ** 2 + deltay ** 2 + deltaz ** 2)) ** 2 / w1 / w1)
  ax = - 2 / m / w / w * (u0 * x * f + u1 * deltax * exp)
  ay = - 2 / m / w / w * (u0 * y * f + u1 * deltay * exp)
  az = - 2 / m / w / w * (u0 * z * f + u1 * deltaz * exp)
  return [vx, ax, vy, ay, vz, az]


def probability():
    count = 0
    for i in range(n):
        xi, yi, zi, vxi ,vyi, vzi = np.random.normal(0, stx, 1)[0],np.random.normal(0, stx, 1)[0],np.random.normal(0, stx, 1)[0],np.random.normal(0, stp, 1)[0] / m,np.random.normal(0, stp, 1)[0] / m,np.random.normal(0, stp, 1)[0] / m
        while m * ( vxi ** 2 + vyi ** 2 + vzi ** 2) / 2 - u0 * np.exp(-(xi ** 2 + yi ** 2) / w / w - zi ** 2 / zr / zr) > 0:
            xi, yi, zi, vxi ,vyi, vzi = np.random.normal(0, stx, 1)[0],np.random.normal(0, stx, 1)[0],np.random.normal(0, stx, 1)[0],np.random.normal(0, stp, 1)[0] / m,np.random.normal(0, stp, 1)[0] / m,np.random.normal(0, stp, 1)[0] / m
        S0 = (xi, vxi, yi, vyi ,zi, vzi)
        sol = spi.odeint(dsdt, y0=S0, t=t, tfirst=True, )
        if (m / 2 * (sol[-1][1] ** 2 + sol[-1][3] ** 2 + sol[-1][5] ** 2) - u0 * np.exp(-(sol[-1][0] ** 2 + sol[-1][2] ** 2) / w / w - sol[-1][4] ** 2 / zr / zr)) < 0:
            count += 1
        #print(i)
    return( count / n)

Y, P_Y = [0.01, 1.5], []

for i in Y:

  alpha, beta = 30, 60
  x_dir, y_dir, z_dir = 1, 1, 1
  R, l,u, m, w, zr = 4, 5, 0.5, 1.45, 1, 1
  u1, u0, kt = 1000, 100, 69

  w1 = i
  print(i, end = ', ')
  k = 0.6 * u0
  stx, stp = (kt / 2 / k) ** 0.5, (m * kt) ** 0.5
  tau, count = 2 * l / u, 0
  t = np.arange(0, tau, dt)

  var = np.sqrt(x_dir **2 + y_dir **2 + z_dir ** 2)
  cosa, cosb, cosg = x_dir / var, y_dir / var, z_dir / var
  x1, y1, z1 = R * np.cos(beta / 180 * 3.14) * np.cos(alpha / 180 * 3.14), R * np.cos(beta / 180 * 3.14) * np.sin(alpha / 180 * 3.14), R * np.sin(beta / 180 * 3.14)
  prob = probability()
  P_Y.append(prob)
for i in P_Y:
  print(i, end = ',')
plt.grid()
plt.title("Зависимость вероятности финитного движения \n от глубины стационарной ямы")
plt.xlabel("u1, 10^(-29) Дж")
plt.ylabel("P(u1)")
plt.plot(Y, P_Y)