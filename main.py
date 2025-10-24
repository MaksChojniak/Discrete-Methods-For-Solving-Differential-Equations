import numpy as np
import sympy as sp
from scipy.integrate import odeint

import matplotlib.pyplot as plt

def F(x, y):
    return (1.3*x) + (1.4*y)


def prosta_metoda_eulera(x0, y0, n, h):
    x = np.zeros(n+1, dtype=float)
    for i in range(n+1):
        x[i] = x0 + (h*i)
    y = np.zeros(n+1, dtype=float)
    y[0] = y0
    for i in range(0, n):
        y[i+1] = y[i] + (h*F(x[i], y[i]))

    return x, y

def ulepszona_metoda_eulera(x0, y0, n, h):
    x = np.zeros(n+1, dtype=float)
    for i in range(n+1):
        x[i] = x0 + (h*i)
    y = np.zeros(n+1, dtype=float)
    y[0] = y0

    for i in range(0, n):
        _x = x[i] + h
        _y = y[i] + (h*F(x[i],y[i]))
        _m = F(_x, _y)
        y[i+1] = y[i] + ((h * (F(x[i], y[i]) + _m)) / 2)

    return x, y

def metoda_rungego_kutty(x0, y0, n, h):
    x = np.zeros(n+1, dtype=float)
    for i in range(n+1):
        x[i] = x0 + (h*i)
    y = np.zeros(n+1, dtype=float)
    y[0] = y0

    for i in range(n):
        m1 = F(x[i], y[i])
        m2 = F(x[i]+(h/2), y[i]+(h*m1/2))
        m3 = F(x[i]+h, y[i]-(h*m1)+(2*h*m2))
        y[i+1] = y[i] + ( h * (m1+(4*m2)+m3) / 6 )
    return x, y

x0 = 0.1
y0 = 0.0
h = 0.2
n = 5

x_p_euler, y_p_euler = prosta_metoda_eulera(x0, y0, n, h)
x_u_euler, y_u_euler = ulepszona_metoda_eulera(x0, y0, n, h)
x_rungegeo_kutty, y_rungegeo_kutty = metoda_rungego_kutty(x0, y0, n, h)


def model_odeint(y, x):
    return (1.4*y) + (1.3*x)

x_odeint = np.linspace(x0, x0+(h*n), n+1)
y_odeint = odeint(model_odeint, y0, x_odeint).reshape(n+1)

print(f'metoda prosta Eulera: \t\t{np.array2string(y_p_euler, formatter={'float_kind': lambda _x: "%.8f" % _x})}')
print(f'metoda Eulera-Cauchy\'ego: \t{np.array2string(y_u_euler, formatter={'float_kind': lambda _x: "%.8f" % _x})}')
print(f'metoda Rungego-Kutty: \t\t{np.array2string(y_rungegeo_kutty, formatter={'float_kind': lambda _x: "%.8f" % _x})}')
print(f'odeint: \t\t\t{np.array2string(y_odeint, formatter={'float_kind': lambda _x: "%.8f" % _x})}')


x = np.linspace(x0, x0+(h*n), 500)
y = odeint(model_odeint, y0, x)

plt.figure(figsize=(8,8))
plt.plot(x, y, label="odeint")
plt.scatter(x_p_euler, y_p_euler, label='metoda prosta Eulera')
plt.scatter(x_u_euler, y_u_euler, label='metoda Eulera-Cauchy\'ego')
plt.scatter(x_rungegeo_kutty, y_rungegeo_kutty, label='metoda Rungego-Kutty')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Porownanie metod rozwiazywania rownan rozniczkowych zwyczajnych")
plt.legend()
plt.grid()
plt.show()



# x = sp.Symbol('x')
# y = sp.Function('y')
# ode = sp.Eq(y(x).diff(x), (1.3*x) + (1.4*y(x)))
# print("Równanie różniczkowe:")
# sol_ics = sp.dsolve(ode, y(x), ics={y(x0): y0})
# rhs_sol = sol_ics.rhs
# f_sol_analitycznie = sp.lambdify(x, rhs_sol, 'numpy')
# x_points = np.linspace(x0, x0+(h*n), 500)
# y_points = f_sol_analitycznie(x_points)
# plt.figure(figsize=(8,8))
# plt.plot(x_points, y_points, label='rozw. analityczne')
# plt.scatter(x_p_euler, y_p_euler, label='metoda prosta Eulera')
# plt.scatter(x_u_euler, y_u_euler, label='metoda Eulera-Cauchy\'ego')
# plt.scatter(x_rungegeo_kutty, y_rungegeo_kutty, label='metoda Rungego-Kutty')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title("Porownanie metod rozwiazywania rownan rozniczkowych zwyczajnych")
# plt.legend()
# plt.grid()
# plt.show()