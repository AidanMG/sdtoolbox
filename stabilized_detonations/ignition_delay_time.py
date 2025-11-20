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
P1 = 500000 
P1_range=np.logspace(5,6,15)  # Initial pressures to consider
T1 = 850
T1_range = np.linspace(1000,1200,15)  # Initial temperatures to consider

T1_mesh, P1_mesh = np.meshgrid(T1_range, P1_range)
q = 'H2:2 O2:1'
mech = 'cti/Lietal_2003.yaml'

gas1 = ct.Solution(mech)
# tend = 1
# # Set up gas object
# gas1 = ct.Solution(mech)
# gas1.TPX = T1,P1,q
# CPout = cpsolve(gas1,t_end=tend,max_step=1e-1)

def cp_ind_time(T1, P1, q, mech, tend=1):
    gas1.TPX = T1,P1,q
    CPout = cpsolve(gas1,t_end=tend,max_step=1e-1)
    return CPout['ind_time']

TI=np.zeros(T1_mesh.shape)
for i in range(T1_mesh.shape[0]):
    for j in range(T1_mesh.shape[1]):
        # Put a limit on the runtime of each iteration
        print('Computing for T1=%.1f K, P1=%.1f Pa...' % (T1_mesh[i,j], P1_mesh[i,j]))
        
        TI[i,j]=cp_ind_time(T1_mesh[i,j], P1_mesh[i,j], q, mech, tend=1)
        print('T1=%.1f K, P1=%.1f Pa, DI=%.4g s' % (T1_mesh[i,j], P1_mesh[i,j], TI[i,j]))

fig, ax = plt.subplots(figsize=(8,6))
subs_series = np.log10(np.linspace(2,10,9,endpoint=True))
norm= LogNorm(vmin=TI.min(), vmax=TI.max())
cs = ax.contour(T1_mesh, P1_mesh/1e6, TI, locator=ticker.LogLocator(subs=range(1,10)), norm=norm, cmap='viridis')

cbar = fig.colorbar(cs, label='Ignition Delay Time (s)')

ax.set_xlabel('Initial Temperature (K)')
ax.set_ylabel('Initial Pressure (MPa)')
ax.set_title('Ignition Delay Time for H2/O2 Mixture')
plt.show()

# CVout = cvsolve(gas1,t_end=tend,max_step=1e-2)
# cv_plot(CVout,major_species='All', maxt=tend)
#cp_plot(CPout,major_species='All', maxt=tend)
# print('Constant-pressure ignition delay time = %.4g s' % CPout['ind_time'])
# print('Constant-volume ignition delay time = %.4g s' % CVout['ind_time'])