from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve
from sdtoolbox.utilities import CJspeed_plot, znd_plot, znd_fileout
import cantera as ct
import numpy as np

from matplotlib import pyplot as plt

temps = np.linspace(200,300,10)
press = np.linspace(15e3,45e3,10) 

keys = ["ind_len_ZND", "ind_time_ZND", "exo_len_ZND", "exo_time_ZND"]

Ts, Ps = np.meshgrid(temps, press)
znd_outputs = {}

def save():
    np.savetxt("Ts.csv", Ts, delimiter=',')
    np.savetxt("Ps.csv", Ps, delimiter=',')
    for k in keys:
        np.savetxt(f"{k}.csv", znd_outputs[k], delimiter=',')

def load():
    Ts = np.loadtxt("Ts.csv", float, delimiter=',')
    Ps = np.loadtxt("Ps.csv", float, delimiter=',')
    for k in keys:
        znd_outputs[k] = np.loadtxt(f"{k}.csv", float, delimiter=',')
    return Ts, Ps, znd_outputs

def plot(X, Y, Z):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    surf = ax.plot_surface(X, Y, Z)
    plt.show()

use_cached = True
if use_cached:
    Ts, Ps, znd_outputs = load()
    plot(Ts, Ps, znd_outputs['ind_len_ZND'])
    plot()
    quit()





for key in keys:
    znd_outputs[key] = np.zeros_like(Ts, dtype=float)
q = 'H2:1 O2:8.5 CH4:4'
mech = 'cti/ffcm1.yaml'


dict_indef = True


# Find CJ speed and related data, make CJ diagnostic plots
for i, P1 in enumerate(press):
    for j, T1 in enumerate(temps):
        print(f"({i}, {j})")
        cj_speed,R2,plot_data = CJspeed(P1,T1,q,mech,fullOutput=True)
        
        # Set up gas object
        gas1 = ct.Solution(mech)
        gas1.TPX = T1,P1,q

        # Find post shock state for given speed
        gas = PostShock_fr(cj_speed, P1, T1, q, mech)

        # Solve ZND ODEs, make ZND plots
        t_end = 1e-6
        for nt in range(100):
            znd_out = zndsolve(gas,gas1,cj_speed,t_end=1e-4,advanced_output=True)
                
            if znd_out['exo_time_ZND'] == 0:
                continue

            for k in keys:
                znd_outputs[k][i][j] = znd_out[k]
            break

save()


    


        

