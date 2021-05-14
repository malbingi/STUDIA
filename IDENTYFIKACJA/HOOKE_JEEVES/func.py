from math import pow, cos, pi
from numba import jit

@jit(nopython=True)
def fun_1(x): #1,1
    f = 2.5*(x[0]**2 - x[1])**2 + (1-x[1])**2
    return f

@jit(nopython=True)
def Rosenbrock(x): #1,1
    f = 100 * (x[1]-x[0]**2)**2 + (1-x[0])**2
    return f

@jit(nopython=True)
def Rastring(x): #0,0 trudne
    f = 20 + x[0]**2 + x[1]**2 - 10 * (cos(2*pi*x[0]) + cos(2*pi*x[1]))
    return f

def call(x, fun):
    if fun == 0:
        return fun_1(x)
    elif fun == 1:
        return Rosenbrock(x)
    else: 
        return Rastring(x)