import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sympy import symbols, diff,cos, sin, lambdify
import sympy as sym
"""
We consider the following problem:
dT/dt + a * dT/dx = d/dx (k dt/dx) + S(x) in [0, 1]
"""

a = 0
k = 1


def back_euler(k, a):
    x = symbols('x')
    T = 200 * sym.exp(10 * (x - 1) ** 2) / sym.exp(10) + 300
    dTdt = k * diff(diff(T,x),x) - a * (diff(T,x))
    dTdx = diff(T,x)
    ddTdx = - diff(dTdx)

    dTdt =  k * ddTdx + ddTdx - a * dTdx

    return lambdify(x, dTdt), lambdify(x, T), lambdify(x, dTdx), lambdify(x, ddTdx) 





def fem_steady_diffusion(N):
    x = np.linspace(0, 1, N)
    h = x[1] - x[0]

    K = np.zeros((N,N))
    F = np.zeros(N)

    [_,_,_,S] = back_euler(k, a)
    # Assembly

    for i in range(N-1):
        Ke = (1/h) * np.array([[1, -1], [-1, 1]])
        xm = 0.5 * (x[i] + x[i + 1])
        Fe = S(xm) * h * 0.5 * np.array([1,1])

        K[i:i+2, i:i+2] += Ke
        F[i:i+2] += Fe
    
    # Dirichlet BC ar x = 0
    K[0,:] = 0
    K[0,0] = 1
    F[0] = 500

    # Neumann BC at x = 1 -> nothing to add (natural)

    T = np.linalg.solve(K, F)
    return x, T

def fem_backward_euler(N, dt, tfinal):
    x = np.linspace(0, 1, N)
    h = x[1] - x[0]

    M = np.zeros((N,N))
    K = np.zeros((N,N))
    F = np.zeros(N)

    # Gauss points

    gauss_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]
    gauss_wts = [1,1]

    # Assembly
    for i in range(N-1):
        Me = h/6 * np.array([[2,1], [1,2]])
        Ke = (1/h) * np.array([[1, -1], [-1, 1]])

        M[i:i+2, i:i+2] += Me
        K[i:i+2, i:i+2] += Ke

        Fe = np.zeros(2)
        [dTTdt, T, dTdx, S] = back_euler(k,a)
        for gp, w in zip(gauss_pts, gauss_wts):
            xg = 0.5 * (x[i + 1] + x[i]) + 0.5 * h * gp
            phi = np.array([0.5 * (1 - gp), 0.5 * (1 + gp)])
            Fe += S(xg) * phi * w * 0.5 * h

        F[i:i+2] += Fe

    # Initial condition
    Tn = np.ones(N) * 300.0

    # Dirichlet BC at x=0
    def apply_dirichlet(A, b):
        A[0,:] = 0.0
        A[0,0] = 1.0
        b[0] = 500.0

    nt = int(tfinal/dt)
    A = M + dt*K

    for _ in range(nt):
        rhs = M @ Tn + dt*F
        apply_dirichlet(A, rhs)
        Tn = np.linalg.solve(A, rhs)

    return x, Tn

def main():

    [dTdt, T, dTdx, _] = back_euler(k, a)
    N = 100

    x , T_num = fem_steady_diffusion(N)
    T_ref =  T(x)

    plt.plot(x, T_num, label = "FEM")
    plt.plot(x, T_ref, '--', label = "Exact")
    plt.xlabel("x")
    plt.ylabel("T")
    plt.legend()
    plt.grid()
    plt.show()

    # Error

    error = np.linalg.norm(T_num - T_ref, 2) / np.sqrt(N)
    print("L2 error =", error)

    N = 200
    dt = 1e-3
    tfinal = 0.1
    x, T_num = fem_backward_euler(N,dt, tfinal)
    T_ref = T(x)

    # Error 
    dx = x[1] - x[0]
    error = np.sqrt(np.sum((T_num - T_ref) ** 2) * dx)
    print("L2 error =", error)

    # Plot
    plt.plot(x, T_num, label = "Backward Euler FEM")
    plt.plot(x, T_ref, '--', label = "Exact steady")
    plt.legend()
    plt.grid()
    plt.show()


main()