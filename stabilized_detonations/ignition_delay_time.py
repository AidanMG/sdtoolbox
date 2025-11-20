from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.cp import cpsolve, CPSys
from sdtoolbox.cv import cvsolve, CVSys
from sdtoolbox.utilities import CJspeed_plot, cp_plot, cv_plot
import cantera as ct
import datetime


P1 = 500000 
T1 = 975
q = 'H2:2 O2:1'
mech = 'cti/Lietal_2003.yaml'

tend = 0.01
# Set up gas object
gas1 = ct.Solution(mech)
gas1.TPX = T1,P1,q
CPout = cpsolve(gas1,t_end=tend,max_step=1e-5)
# CVout = cvsolve(gas1,t_end=tend,max_step=1e-2)
# cv_plot(CVout,major_species='All', maxt=tend)
cp_plot(CPout,major_species='All', maxt=tend)
print('Constant-pressure ignition delay time = %.4g s' % CPout['ind_time'])
# print('Constant-volume ignition delay time = %.4g s' % CVout['ind_time'])