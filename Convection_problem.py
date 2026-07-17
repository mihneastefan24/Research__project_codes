#%%
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import os

os.chdir("C:/Users/mihne/Desktop/scoala/VKI/VKI/NLAB")
# %%
plt.close("all")

# -----------------------------
# Parameters
# -----------------------------
L = 1.0
kappa = 225.0          # thermal conductivity
tfinal = 1e-3
T_left = 800.0         # Dirichlet BC at x=0
# Nx_vector = np.array([41, 51, 101, 201, 251])
Nx_vector = np.array([41, 51, 101, 201, 251])


def laplacian(T, dx, kappa):
    d2T = np.zeros_like(T)
    N = len(T)

    # ---- Enforce Neumann BC at right boundary (4th-order) ----
#    T[-1] = ( 48*T[-2] - 36*T[-3] + 16*T[-4] - 3*T[-5]) / 25.0
#    T[-1] = (- T[-3] + 16 * T[-2] - 30 * T[-1] + 16 * T[0] - T[1]) / (12 * dx ** 2)
    # ---- Interior (central 4th-order) ----
    d2T[2:-2] = (
        -T[4:] + 16*T[3:-1] - 30*T[2:-2]
        + 16*T[1:-3] - T[0:-4]
    ) / (12 * dx**2)

# Forth order in the boundary layer leads to instabilities

    # ---- Left boundary closure (i=1) ----
#    d2T[1] = (        -T[3]        + 16*T[2]        - 30*T[1]        + 16*T[0]        - 1*T[-1]    ) / (12 * dx**2)
    d2T[1] = ( T[2] - 2* T[1] + T[0] ) / dx ** 2
    # ---- Right boundary closure (i=N-2) ----
#    d2T[-2] = (        -T[0]        + 16*T[-1]        - 30*T[-2]        + 16*T[-3]        - 1*T[-4]   ) / (12 * dx**2)
    d2T[-2] = (T[-1] - 2 * T[-2] + T[-3]) / dx ** 2
    # Dirichlet node (second derivative not used)
    d2T[0] = 0.0
    d2T[-1] = d2T[-2]

    return kappa * d2T

def heat_rhs(T, dx, source):

    d2T = np.zeros_like(T)
    d2T[1:-1] = (-T[1:] + 16 * T[2:0] - 30 * T[1:-1] + 16 * T[0:-2] - T[-1:-3]) / (dx ** 2 * 12)
    d2T[-1] = d2T[-2]
    return kappa * d2T + source

def solve_physical_heat(Nx, tf, kappa):
    dx = L / (Nx - 1)
    dt = 1/24 * dx ** 2 / kappa
    Nt = int(tf / dt)

    x = np.linspace(0, L, Nx)
    # Initial condition
    T = np.ones(Nx) * 300.0
    T[0] = T_left

    # Time integration

    for _ in range(Nt):
        T[0] = T_left
        T[-1] = T[-2]

        k1 = laplacian(T,dx, kappa)
        T_star = T + dt * k1

        T_star[0] = T_left
        T_star[-1] = T_star[-2]

        k2 = laplacian(T_star, dx, kappa)
        T_star = T_star + dt * k2

        T_star[0] = T_left
        T_star[-1] = T_star[-2]

        k3 = laplacian(T_star,dx, kappa)
        T_star = T_star + dt * k3

        T_star[0] = T_left
        T_star[-1] = T_star[-2]

        k4 = laplacian(T_star, dx, kappa)
        T += 1/6 * dt * (k1 + 2 * k2 + 2 * k3 + k4)

        # Enforce BCs
        T[0] = T_left
        T[-1] = T[-2]

    return x, T


plt.figure(1)

for Nx in Nx_vector:

    x, T = solve_physical_heat(Nx, tfinal, kappa)

    plt.plot(x, T, label = f"dx is {L/(Nx-1):.4f}")
    plt.xlabel("x")
    plt.ylabel("Temperature [K]")
    plt.title("1D Heat Equation Solution")
    plt.grid()
    
plt.legend()
plt.show() 

# %% 
## ==================================
#  Create the manufactured solution
## ==================================

tfinal = 1e-3
def periodic_laplacian(T,dx):
    im2 = np.roll(T, 2)
    im1 = np.roll(T, 1)
    ip1 = np.roll(T, -1)
    ip2 = np.roll(T, -2)
    return (-ip2 + 16*ip1 - 30*T + 16*im1 - im2) / (12 * dx**2)

def solve_heat(Nx, T_exact_func, S_func):
    x = np.linspace(0, 1.0, Nx, endpoint = False)
    dx = x[1] - x[0]
    kappa = 225

    # Stable timestep
    nt = int(1e5)
    dt = tfinal / nt
    t = np.linspace(0, tfinal, nt)

    # Initial condition from MMS
    T = T_exact_func(x, 0.0)

    for n in range(nt):
        t = n * dt
        S = S_func(x, t)

        def RHS(U):
            return kappa * periodic_laplacian(U.copy(), dx) + S
            
        k1 = dt * RHS(T)
        k2 = dt * RHS(T + 0.5 * k1)
        k3 = dt * RHS(T + 0.5 * k2)
        k4 = dt * RHS(T + k3)

        T = T + (1/6)*(k1 + 2 * k2 + 2* k3 +  k4)

        # Enforce BCs
        T[0] = 800.0
    return x, T

x_sym, t_sym = sym.symbols('x t')
steady = np.array([True, False]) 
for s_i in range(1):

    if steady[s_i]:
        T_mms = 800 -  sym.sin(2*sym.pi*x_sym/L)
        source_term = sym.diff(T_mms,t_sym) -kappa * sym.diff(T_mms,x_sym,2)  # Version of steady state
    else:
        T_mms = 800 - sym.exp(-kappa * (sym.pi / (2 * L)) ** 2 * t_sym) * sym.sin(sym.pi * x_sym / (2 * L))
        source_term = -kappa * sym.diff(T_mms, x_sym, 2) + sym.diff(T_mms, t_sym, 1)


    S_sym = sym.simplify(source_term)

    print("Manufactured solution :", T_mms)
    print("Source term: ", S_sym)

    T_exact_func = sym.lambdify((x_sym, t_sym), T_mms, "numpy")
    S_func = sym.lambdify((x_sym, t_sym), S_sym, "numpy")

    # =====================================================
    # Convergence Study
    # =====================================================
    Nx_list = [41, 51,161, 201, 251]
    errors = []
    dx_list = []

    for Nx in Nx_list:
        x, T_num = solve_heat(Nx, T_exact_func, S_func)
        dx = x[1] - x[0]

        T_exact = T_exact_func(x, tfinal)
        #decide whetere you continue with the sqrt(Nx) or not
        error = np.linalg.norm(T_num - T_exact, 2) * np.sqrt(dx)

        errors.append(error)
        dx_list.append(dx)

    # Compute observed order
    order = np.polyfit(np.log(dx_list), np.log(errors), 1)
    print("Observed order:", -order[0])

    # Plot convergence
    plt.figure()
    plt.loglog(dx_list, errors, 'o-', label="Error")
    plt.loglog(dx_list,
            errors[0]*(np.array(dx_list)/dx_list[0])**4,
            '--', label= r"O($dx^4$)")
    plt.gca().invert_xaxis()
    plt.legend()
    plt.grid(True)
    plt.xlabel("dx [m]")
    plt.ylabel("L2 Error")
    plt.title("Convergence Study")
    plt.show()

    mesh = np.linspace(0, L, Nx)
    plt.figure()
    plt.plot(mesh, T_num,label=r"$T_{numerical}$",color='r', linewidth = 3.0)
    plt.plot(mesh, np.array(T_exact), label = r"$T_{reference}$", linestyle = ':', linewidth = 3.0)
    plt.legend()
    plt.title("Temperature evolution comparison")
    plt.grid(True)
    plt.xlabel("x [m]")
    plt.ylabel("Temperature [K]")
    #plt.ylim([750, 850])
    plt.show()


#%% Implement a variable k(x)

# =====================================================
# Parameters
# =====================================================

L = 1.0
tfinal = 1e-3
kappa = 225


# =====================================================
# 4th order periodic first derivative
# =====================================================

def D1_4th(T,dx):

    im2 = np.roll(T,2)
    im1 = np.roll(T,1)
    ip1 = np.roll(T,-1)
    ip2 = np.roll(T,-2)

    return (-ip2 + 8*ip1 - 8*im1 + im2)/(12 * dx)

def periodic_laplacian(T,dx):


    im2 = np.roll(T, 2)
    im1 = np.roll(T, 1)
    ip1 = np.roll(T, -1)
    ip2 = np.roll(T, -2)
    return (-ip2 + 16*ip1 - 30*T + 16*im1 - im2) / (12 * dx**2)



def solve_heat_k(Nx, k_func, T_func, S_func):
    x = np.linspace(0, L, Nx, endpoint=False)
    dx = x[1] - x[0]
    dt = 0.1 * dx **2 / 225
#    dt = tfinal / nt
    nt = int(tfinal / dt)
    T = T_func(x)
    kappa = k_func(x)
    
    for n in range(nt):
        RHS = periodic_laplacian(T, dx) * kappa + D1_4th(T, dx) * D1_4th(kappa, dx) + S_func(x)
        T = T + dt * RHS
    
    return x, T

x_sym, t_sym = sym.symbols('x t')
k_sym = 225 * (1 + 0.3 * sym.sin( 2 *sym.pi * x_sym / L))
T_mms = 800 - sym.sin(2 * sym.pi * x_sym / L)
S_sym = sym.diff(T_mms, t_sym, 1) - sym.diff(k_sym * sym.diff(T_mms, x_sym, 1), x_sym, 1)

S_sym = sym.simplify(S_sym)

print("Manufactured solution:")
print(T_mms)

print("\nSource term:")
print(S_sym)

print("\nK(x):")
print(k_sym)

T_func = sym.lambdify((x_sym), T_mms, "numpy")
S_func = sym.lambdify((x_sym), S_sym, "numpy")
k_func = sym.lambdify((x_sym), k_sym, "numpy")

Nx_vector = [21, 51, 101, 151, 201, 251]

errors = []
dx_vector = []


for Nx in Nx_vector:

    x = np.linspace(0, 1, Nx, endpoint= False)
    dx = x[1] - x[0]
    T_exact = T_func(x)
    x, T_num = solve_heat_k(Nx, k_func, T_func,S_func)
    error = np.linalg.norm(T_num - T_exact) / np.sqrt(Nx) 
    errors.append(error)
    dx_vector.append(dx)
    plt.figure()    
    plt.plot(x,T_num,label=f"Numerical for Nx={Nx}",linewidth=3)
    plt.plot(x,T_exact,':',label=r"Exact",linewidth=3)     
    plt.ylabel("T [K]")
    plt.xlabel("x [m]")
    plt.title("Difference ")
    plt.legend()
    plt.grid()
    plt.show()



orderA = np.polyfit(np.log(dx_vector), np.log(errors),1)



#%%
# =====================================================
# Plot convergence
# =====================================================

plt.figure()

plt.loglog(dx_vector,errors,'o-',label="Error")

plt.loglog(
dx_vector,
errors[0]*(np.array(dx_vector)/dx_vector[0])**4,
'--',
label=r"O($dx^4$)"
)

plt.gca().invert_xaxis()

plt.xlabel("dx [m]")
plt.ylabel("L2 Error")
plt.title("Error evolution for variable k")
plt.legend()
plt.grid(True)

plt.show()


# =====================================================
# Final Solution Plot
# =====================================================

plt.figure()

plt.plot(x,T_num,label="Numerical",linewidth=3)
plt.plot(x,T_exact,':',label="Exact",linewidth=3)
plt.ylabel("T [K]")
plt.xlabel("x [m]")
plt.title("Numerical vs Exact evolution for variable k")
plt.legend()
plt.grid()
plt.show()


#%%
# =====================================================
# Time evolution toward steady state
# Temperature should go from 300 K → 800 K
# =====================================================

import numpy as np
import matplotlib.pyplot as plt

L = 1.0
kappa = 225.0
T_left = 800.0
Nx = 101
tfinal = 0.02   # longer time to see convergence

dx = L/(Nx-1)
dt = 0.2*dx**2/kappa
Nt = int(tfinal/dt)

x = np.linspace(0,L,Nx)

# Initial condition
T = np.ones(Nx)*300.0

time_vec = []
T_mid_vec = []

mid_index = Nx//2


def laplacian_3rd(T, dx):

    d2T = np.zeros_like(T)

    d2T[2:-2] = (
        -T[4:] + 16*T[3:-1] - 30*T[2:-2]
        + 16*T[1:-3] - T[0:-4]
    )/(12*dx**2)

    d2T[1] = (T[2]-2*T[1]+T[0])/dx**2
    d2T[-2] = (T[-1]-2*T[-2]+T[-3])/dx**2

    d2T[0]=0
    d2T[-1]=d2T[-2]

    return kappa*d2T


# Time integration (midpoint scheme)

for n in range(Nt):

    t = n*dt

    T[0] = T_left
    T[-1] = T[-2]

    k1 = laplacian_3rd(T,dx)

    Tstar = T + 0.5*dt*k1

    Tstar[0] = T_left
    Tstar[-1] = Tstar[-2]

    k2 = laplacian_3rd(Tstar,dx)

    T = T + dt*k2

    time_vec.append(t)
    T_mid_vec.append(T[mid_index])


# Plot

plt.figure()

plt.plot(time_vec,T_mid_vec,linewidth=3)

plt.xlabel("Time [s]")
plt.ylabel("Temperature at x=L/2 [K]")

plt.title("Time Evolution Toward Steady State")

plt.grid()

plt.show()