from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.cp import cpsolve, CPSys
from sdtoolbox.cv import cvsolve, CVSys
from sdtoolbox.utilities import CJspeed_plot, cp_plot, cv_plot
import cantera as ct
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import os

def md_loc(eta_0):
    return 0.67*np.sqrt(eta_0)

P1_range=np.logspace(4,7,20)  # Initial pressures to consider
T1_range = np.linspace(1000,1400,20)  # Initial temperatures to consider
T1_mesh, P1_mesh = np.meshgrid(T1_range, P1_range)
q = 'H2:2 O2:1'
mech = 'cti/Lietal_2003.yaml'
gas1 = ct.Solution(mech)


def cp_ind_time(T1, P1, q, tend=1):
    gas1.TPX = T1,P1,q
    CPout = cpsolve(gas1,t_end=tend,max_step=1e-1)
    return CPout['ind_time']

def cp_grid(T1_mesh, P1_mesh, q):

    # tend = 1
    # # Set up gas object
    # gas1 = ct.Solution(mech)
    # gas1.TPX = T1,P1,q
    # CPout = cpsolve(gas1,t_end=tend,max_step=1e-1)



    TI=np.zeros(T1_mesh.shape)

    for i in range(T1_mesh.shape[0]):
        for j in range(T1_mesh.shape[1]):
            # Put a limit on the runtime of each iteration
            print('Computing for T1=%.1f K, P1=%.1f Pa...' % (T1_mesh[i,j], P1_mesh[i,j]))
            
            TI[i,j]=cp_ind_time(T1_mesh[i,j], P1_mesh[i,j], q, tend=1)
    return TI


if os.path.exists('ignition_delay_time.npz'):
    data = np.load('ignition_delay_time.npz')
    T1_mesh = data['T1_mesh']
    P1_mesh = data['P1_mesh']
    TI = data['TI']
else:
    TI = cp_grid(T1_mesh, P1_mesh, q)
    np.savez('ignition_delay_time.npz', T1_mesh=T1_mesh, P1_mesh=P1_mesh, TI=TI)
# Using velocity Mach 1 roughly 540 m/s
U1=540
fig, ax = plt.subplots(figsize=(8,6))
subs_series = np.log10(np.linspace(2,10,9,endpoint=True))
norm= LogNorm()
cs = ax.contourf(T1_mesh, P1_mesh, TI*U1-1e-3*md_loc(P1_mesh/101325), cmap='viridis', norm=norm, locator=ticker.LogLocator(subs='all'))
cb = fig.colorbar(cs, ticks = ticker.LogLocator(subs=range(10)))
cb.ax.minorticks_on()
ax.set_yscale('log')
cb.set_label('Ignition Delay length (m)')
ax.set_xlabel('Initial Temperature (K)')
ax.set_ylabel('Initial Pressure (Pa)')
ax.set_title(f'Pre-shock Ignition Delay length for {q} using {mech} mechanism')
plt.show()


# CVout = cvsolve(gas1,t_end=tend,max_step=1e-2)
# cv_plot(CVout,major_species='All', maxt=tend)
#cp_plot(CPout,major_species='All', maxt=tend)
# print('Constant-pressure ignition delay time = %.4g s' % CPout['ind_time'])
# print('Constant-volume ignition delay time = %.4g s' % CVout['ind_time'])