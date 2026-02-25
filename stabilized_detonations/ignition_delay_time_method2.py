import pandas as pd
from plotting_utilities import extract_paraview_csv, get_M, average_last_timesteps
import numpy as np

from centerline_model import get_full_model

from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.cp import cpsolve, CPSys
from sdtoolbox.cv import cvsolve, CVSys
from sdtoolbox.utilities import CJspeed_plot, cp_plot, cv_plot
import cantera as ct
import datetime
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import os
from scipy.optimize import fsolve


def method2(gas, d, x_end, N=1000, method="cv", map=None):
    Nsp = gas.n_species
    
    if method=="cp":
        solver = lambda gas, t_end, max_step: cpsolve(gas, t_end, max_step, skip_eval=True)
    elif method=="cv":
        solver = lambda gas, t_end, max_step: cvsolve(gas, t_end, max_step, skip_eval=True)
    else:
        raise Exception(f"Unrecognized solver method: {method} in chemistry step")
    
    
    
    
    xs = np.linspace(0,1,N)
    if map is not None:
        xs = map(xs)
    xs = xs * x_end
        
    
    M_func = get_full_model() # M(x/d)
    
    Ms = M_func(xs/d)
    Ts = np.zeros_like(xs)
    Ps = np.zeros_like(xs)
    Ys = np.zeros(shape=(N,Nsp))
    Hs = np.zeros_like(xs)
    Ss = np.zeros_like(xs)
    ts = np.zeros_like(xs)
    Mols = np.zeros_like(Ys)
    dTdt = np.zeros_like(Ys)
    
    
    ts[0] = 0
    Ts[0] = gas.T
    Ps[0] = gas.P
    Ys[0] = gas.Y
    Mols[0] = gas.X
    Hs[0] = gas.enthalpy_mass + 0.5*(Ms[0] * gas.sound_speed)**2 
    Ss[0] = gas.entropy_mass
    dTdt[0] = 0
    
    
    for i in range(1,len(xs)):
        gas = isentropic(gas, Ms[i-1], Ms[i])   # Update gas properties to new state
        ux = Ms[i] * gas.sound_speed              # Current cell velocity
        dx = (xs[i]-xs[i-1])
        
        dt = dx / ux # u = dx/dt , dt = dx / u
        
        out = solver(gas, dt, dt/100)
        if method == "cp":
            Ps[i] = gas.P
            Ts[i] = out['T'][-1]
        else:
            Ts[i], Ps[i] = out['T'][-1], out['P'][-1]
        Ys[i] = out['speciesY'].T[-1]
        Mols[i] = out['speciesX'].T[-1]
        ts[i] = ts[i-1] + dt
        gas.TPY = Ts[i], Ps[i], Ys[i]
        Hs[i] = gas.enthalpy_mass + 0.5*(Ms[i] * gas.sound_speed)**2 
        Ss[i] = gas.entropy_mass
        dTdt[i] = out['dTdt'][-1]
        
    solution = {
        "x": xs,
        "t": ts,
        "M": Ms,
        "T": Ts,
        "P": Ps,
        "Y": Ys,
        "X": Mols,
        "h_stag": Hs,          # total (stagnation) enthalpy
        "s": Ss,
        "dTdt": dTdt,          # Chemical dTdt only
        "species_names": gas.species_names
    }
    
    return solution

def isentropic(gas, M0, M1):
    """Expands/compresses a gas object isentropically from state M=M0 to M=M1

    Finds Temperature and Pressure of M1 state by enforcing:
        deltaS = (s(M0) - s(M1)) = 0                    -> Entropy Constraint 
        delta(H_stag) = (h_stag(M0)-h_stag(M1)) = 0     -> Energy Conservation
        
    
    Args:
        gas (ct.Solution): gas object at TP -> M0
        M0 (float): M0 
        M1 (float): M1
    """
    # ---- Initial State ---- 
    h0 = gas.enthalpy_mass
    s0 = gas.entropy_mass
    a0 = gas.sound_speed
    
    # Calculate the stagnation enthalpy
    h_stag = h0 + 0.5 * (M0 * a0)**2 
    
    def residual_TP(X):
        T, P = X
        gas.TP = T, P
        s = gas.entropy_mass
        h = gas.enthalpy_mass
        a = gas.sound_speed
        return [
            s - s0,
            h + 0.5*(M1*a)**2 - h_stag
        ]

    T, P = fsolve(residual_TP, [gas.T, gas.P], epsfcn = 1e-8)
    gas.TP = T, P
    
    return gas
    
def plot_species_vs_time(solution, species=None, logy=False):
    """
    Plot species mass fractions Yi vs time.

    Parameters
    ----------
    solution : dict
        Output dictionary from method2().
        Must contain keys: 't', 'Y', 'species_names'.
    species : list of str or None
        Species names to plot. If None, plot all species.
    logy : bool
        If True, use logarithmic y-axis.
    """


    t = solution["t"]
    Y = solution["Y"]
    names = solution["species_names"]

    if species is None:
        species = names

    fig, ax = plt.subplots()

    for sp in species:
        if sp not in names:
            raise ValueError(f"Species '{sp}' not found in solution.")
        idx = names.index(sp)
        ax.plot(t, Y[:, idx], label=sp)

    if logy:
        ax.set_yscale("log")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mass Fraction")
    ax.grid(True)

    # Place legend outside
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    
    # Make room for legend
    fig.tight_layout()
    fig.subplots_adjust(right=0.75)

    plt.show()
    


def expansion_model(t, dia):
    
    pass
    

if __name__ == "__main__":
    pass
    gas = ct.Solution("cti/Lietal_2003.yaml")
    gas.TPX = 1300, 2*1013250, "H2:2, O2:1"
    dia = 0.025
    xmax = dia*5
    
    power_law_map = lambda xi, alpha: xi**alpha
    
    sol = method2(gas, dia, xmax, N=100, method='cv', map=lambda x: power_law_map(x, alpha=4))
    # plt.plot(sol['t'], sol['dTdt'])
    plt.plot(sol['t'], sol['P'])
    plt.show()
    plot_species_vs_time(sol, logy=False)
    
    
    
    
    # #mach_profile = get_full_model()
    # T0, P0 = 1000, 1013250
    # q = "H2:2 O2:1"

    # gas = ct.Solution("cti/Lietal_2003.yaml")

    # # -------------------------
    # # 1. Stagnation state
    # # -------------------------
    # gas.TPX = T0, P0, q

    # h0 = gas.enthalpy_mass
    # s0 = gas.entropy_mass

    # # -------------------------
    # # 2. Define residual
    # # -------------------------

    # def residual(T):

    #     # For given T, solve for P from isentropy
    #     def s_residual(P):
    #         gas.TP = T, P
    #         return gas.entropy_mass - s0

    #     sol = root_scalar(s_residual, bracket=[1e2, P0*10])
    #     P = sol.root

    #     gas.TP = T, P

    #     h = gas.enthalpy_mass
    #     a = gas.sound_speed

    #     return h + 0.5*a**2 - h0

    # # -------------------------
    # # 3. Solve for T
    # # -------------------------
    
    # solT = root_scalar(residual, bracket=[T0*0.1, T0])
    # T_static = solT.root

    # # Compute final pressure
    # def s_residual(P):
    #     gas.TP = T_static, P
    #     return gas.entropy_mass - s0

    # solP = root_scalar(s_residual, bracket=[1e2, P0*10])
    # P_static = solP.root

    # print(T_static, P_static)
    
    
    
    

