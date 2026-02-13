import pandas as pd
from plotting_utilities import extract_paraview_csv, get_M, average_last_timesteps
import numpy as np


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



def method1(q, csv_filepath="PR20.csv", nT=1000):
    q = ""
    

if __name__ == "__main__":
    # Get data for PR20 run
    data20 = extract_paraview_csv("PR20.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length", "Gamma", "Temperature"])
    data20 = get_M(data20)
    
    Dia = 0.1
    
    x, y, _ = average_last_timesteps(data20, "arc_length", "M")
    xG, yG, _= average_last_timesteps(data20, "arc_length", "Gamma")
    xT, yT, _ = average_last_timesteps(data20, "arc_length", "Temperature")
    T0=yT[0]*(1+(yG[0]-1)/2*(y[0])**2)
    T_s = T0/(1+((yG-1)/2)*y**2)
    
    
    
    
    
    

        
    pass

