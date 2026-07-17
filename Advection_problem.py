#%%
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

plt.close("all")

# =====================================================
# Problem parameters
# =====================================================

L = 1.0           # domain length
k = 1             # conductivity
a = 0             # convection coefficient
tfinal = 1       # final time
T0 = 500          # Dirichlet BC at x = 0

# =====================================================
# Exact analytical solution (manufactured solution)
# =====================================================

x_sym = sym.symbols('x')

T_sym = 200 * sym.exp(10*(x_sym-1)**2)/np.exp(10) + 300
T_d_sym = sym.diff(T_sym, x_sym, 1)
T_dd_sym = sym.diff(T_sym, x_sym, 2)

# Source term obtained from PDE
S_sym = a*T_d_sym - k*T_dd_sym

# Convert symbolic functions to numerical ones
T_func = sym.lambdify(x_sym, T_sym, 'numpy')
T_d_fun = sym.lambdify(x_sym, T_d_sym, 'numpy')
S_func = sym.lambdify(x_sym, S_sym, 'numpy')

# =====================================================
# PART 1 : Steady FEM solution
# =====================================================

N = 10                       # number of elements
x = np.linspace(0, L, N+1)   # mesh
h = L/N

# FEM matrices
K = np.zeros((N+1,N+1))
F = np.zeros(N+1)

# ---- Element loop ----
for e in range(N):

    x1 = x[e]
    x2 = x[e+1]
    he = x2 - x1

    # Local stiffness matrix
    Ke = (k/he)*np.array([[1,-1],
                          [-1,1]])

    # Load vector
    Fe = np.zeros(2)

    # Gauss points
    xi = [-1/np.sqrt(3), 1/np.sqrt(3)]

    for g in xi:

        # mapping to physical element
        xg = 0.5*(x1+x2) + 0.5*he*g

        # shape functions
        N1 = (x2-xg)/he
        N2 = (xg-x1)/he

        Fe += S_func(xg)*np.array([N1,N2])

    Fe *= he/2

    # Assembly
    K[e:e+2,e:e+2] += Ke
    F[e:e+2] += Fe

# ---- Dirichlet BC ----
K[0,:] = 0
K[0,0] = 1
F[0] = T0

# Solve linear system
T = np.linalg.solve(K,F)

# Exact solution
T_exact = T_func(np.linspace(0,1,1000, endpoint = True))

# ---- Plot steady solution ----
plt.figure()
plt.plot(x,T,'o-',label='FEM')
plt.plot(np.linspace(0,1,1000,endpoint = True),T_exact,label='Exact')
plt.xlabel("x[m]")
plt.ylabel("Temperature [K]")
plt.legend()
plt.grid()
plt.show()


#%%
# ---- Error ----
error = np.linalg.norm(T-T_func(x))/np.sqrt(len(x))
print("L2 error =",error)

print("====================")
print("The source term is:")
print(S_sym)
print("====================")
print("Vector Form: ")
print(S_func(x))

#%%
# =====================================================
# PART 2 : Transient problem (Backward Euler)
# =====================================================

# Mass matrix
M = np.zeros((N+1,N+1))

# Rebuild matrices
K = np.zeros((N+1,N+1))
F = np.zeros(N+1)

for e in range(N):

    x1 = x[e]
    x2 = x[e+1]
    he = x2-x1

    Ke = (k/he)*np.array([[1,-1],
                          [-1,1]])

    Me = (he/6)*np.array([[2,1],
                          [1,2]])

    Fe = np.zeros(2)

    xi = [-1/np.sqrt(3),1/np.sqrt(3)]

    for g in xi:

        xg = 0.5*(x1+x2)+0.5*he*g

        N1 = (x2-xg)/he
        N2 = (xg-x1)/he

        Fe += S_func(xg)*np.array([N1,N2])

    Fe *= he/2

    K[e:e+2,e:e+2] += Ke
    M[e:e+2,e:e+2] += Me
    F[e:e+2] += Fe

# ---- Initial condition ----
T = np.ones(N+1)*300
picture_time = [0.01, 0.05, 0.1, 0.5, 1]
count = 0
# Time step
dt = 0.01
time = 0

A = M + dt*K
plt.figure()
marker_vec = ['o','d','^','.','1']
count_marker = 0
while time < tfinal:

    rhs = M@T + dt*F

    # Dirichlet BC
    A[0,:] = 0
    A[0,0] = 1
    rhs[0] = T0

    T = np.linalg.solve(A,rhs)

    time += dt
    if time >= picture_time[count]:
        plt.plot(x,T,'-',marker = marker_vec[count_marker],label=f"Transient at t={np.round(time,2)}")
        count += 1
        count_marker += 1

# Plot transient vs steady
plt.plot(x,T_func(x),label="Steady")
plt.xlabel('x[m]')
plt.ylabel('Temperature[K]')
plt.legend()
plt.grid()
plt.show()

error = np.linalg.norm((T - T_func(x))) / np.linalg.norm(T)

#%%
#%%
# =====================================================
# PART 3 : Time step convergence study
# =====================================================

time_steps = [0.1, 0.05, 0.01, 0.005]
errors = []
count_marker = 0
plt.figure()

for dt in time_steps:

    T = np.ones(N+1)*300
    time = 0
    A = M + dt*K

    while time < tfinal:

        rhs = M @ T + dt * F

        A[0,:] = 0
        A[0,0] = 1
        rhs[0] = T0

        T = np.linalg.solve(A, rhs)

        time += dt

    # compute error
    error = np.linalg.norm(T - T_func(x)) / np.sqrt(len(x))
    errors.append(error)

    # plot temperature profile for this dt
    plt.plot(x, T, marker = marker_vec[count_marker], label=f"Δt = {dt}")
    count_marker += 1

# exact solution
plt.plot(x, T_func(x), 'k--', label="Exact solution")

plt.xlabel("x")
plt.ylabel("Temperature")
#plt.title("Temperature distribution for different time steps")
plt.legend()
plt.grid()
plt.show()

print("Δt =", time_steps)
print("L2 errors =", errors)

# error convergence plot
plt.figure()
plt.loglog(time_steps, errors, 'o-')
plt.xlabel("Δt")
plt.ylabel("L2 error")
plt.grid()
plt.show()
#%%
# =====================================================
# PART 4 : Adaptive mesh refinement
# =====================================================

beta = 3          # mesh grading factor
N = 10            # number of elements

xi = np.linspace(0,1,N+1)
x = xi**beta      # refined mesh near x=1

K = np.zeros((N+1,N+1))
M = np.zeros((N+1,N+1))
F = np.zeros(N+1)

for e in range(N):

    x1 = x[e]
    x2 = x[e+1]
    he = x2-x1

    Ke = (k/he)*np.array([[1,-1],
                          [-1,1]])

    Me = (he/6)*np.array([[2,1],
                          [1,2]])

    Fe = np.zeros(2)

    gp = [-1/np.sqrt(3),1/np.sqrt(3)]

    for g in gp:

        xg = (x1+x2)/2 + he/2*g

        N1 = (x2-xg)/he
        N2 = (xg-x1)/he

        Fe += S_func(xg)*np.array([N1,N2])

    Fe *= he/2

    K[e:e+2,e:e+2] += Ke
    M[e:e+2,e:e+2] += Me
    F[e:e+2] += Fe

# transient solve on refined mesh
T = np.ones(N+1)*300
dt = 0.01
time = 0

A = M + dt*K

while time < tfinal:

    rhs = M@T + dt*F

    A[0,:] = 0
    A[0,0] = 1
    rhs[0] = T0

    T = np.linalg.solve(A,rhs)

    time += dt

# ---- Plot refined mesh result ----
x_exact = np.linspace(0,1,200)

plt.figure()
plt.plot(x,T,'o-',label="Adaptive FEM")
plt.plot(x_exact,T_func(x_exact),label="Exact")
plt.xlabel("x[m]")
plt.ylabel("Temperature[K]")
plt.legend()
plt.grid()
plt.show()