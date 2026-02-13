from scipy.optimize import curve_fit
from plotting_utilities import extract_paraview_csv, get_M, average_last_timesteps
from matplotlib import pyplot as plt
import numpy as np

def nearfield_model(x, y):
    # f = lambda x, A, B: 1 + A*x**B
    f = lambda x, A, B, C: 1 + A*x + B*x**2 + C*x**3
    res = curve_fit(f, x, y, [0.5, 0.5, 0.5], maxfev=10000)
    # print(res[0])
    fret = lambda x: f(x, *res[0])
    return fret

def mixed_model(x, y):
    f = lambda xp, A, B: A*nearfield_model(x,y)(xp) + B*farfield_model_var(x,y)(xp)
    res = curve_fit(f, x, y, [0.5, 0.5])
    fret = lambda x: f(x, *res[0])
    
    f2 = lambda x, A, B, C, D: A + B*x + C*x**2 + D*x**3
    res = curve_fit(f2, x, y, [0.5, 0.5, 0.5, 0.5])
    print(res[0])

    fret = lambda x: f2(x, *res[0])
    return fret
    

def farfield_model(x, y):
    # FIXED G=Gamma MODEL DOES NOT WORK AT ALL
    G=1.4
    f = lambda x, x0, A, C: A*(x-x0)**(G-1) - 0.5*(G+1)/(G-1)/(A*(x-x0)**(G-1))  +  C*(x-x0)**(-3*(G-1)) 
    res = curve_fit(f, x, y, [0.5, 0.5, 0.5], maxfev=10000)
    
    # FLOATING Gamma MODEL
    # f = lambda x, x0, A, C, G: A*(x-x0)**(G-1) - 0.5*(G+1)/(G-1)/(A*(x-x0)**(G-1)) # + C*(x-x0)**(-3*(G-1))
    # res = curve_fit(f, x, y, [-0.5, 0.5, 0.5, 1.4], maxfev=10000)    
    fret = lambda x: f(x, *res[0])
    return fret


def farfield_model_var(x, y):
    G=1.4
    x0 , A, C = 0.40, 3.65, 0.20 
    # f = lambda x, k1, k2, k3: (A*(x-x0)**(G-1))**k1 - (0.5*(G+1)/(G-1)/(A*(x-x0)**((G-1))))**k2 + (C*(x-x0)**(-3*(G-1)))**k3
    f = lambda x, k1, k2: k1*(A*(x-x0)**(G-1) - 0.5*(G+1)/(G-1)/(A*(x-x0)**(G-1))  +  C*(x-x0)**(-3*(G-1)))**k2
    res = curve_fit(f, x, y, [0.5, 0.5], maxfev=10000)
    print("AS Coeffs:", res[0])
    fret = lambda x: f(x, *res[0])
    return fret



if __name__ == "__main__":
    data20 = extract_paraview_csv("PR20.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length", "Gamma", "Temperature"])
    smD20 = extract_paraview_csv("smD_PR20.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length", "Gamma", "Temperature"])
    smD20 = get_M(smD20)
    data20 = get_M(data20)
    bigD20 = extract_paraview_csv("bigD_PR20.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length", "Gamma", "Temperature"])
    bigD20 = get_M(bigD20)    
    
    smDia = 0.01
    Dia = 0.1
    bigDia = 1
    
    smx, smy, _ = average_last_timesteps(smD20, "arc_length", "M", n_last=1)
    bigx, bigy, _ = average_last_timesteps(bigD20, "arc_length", "M", n_last=1)
    
    x, y, _ =average_last_timesteps(data20, "arc_length", "M")
    xG, yG, _= average_last_timesteps(data20, "arc_length", "Gamma")
    xT, yT, _ = average_last_timesteps(data20, "arc_length", "Temperature")

    x = x/Dia # Normalize wrt inlet dia
    smx = smx/smDia
    bigx = bigx/bigDia
    
    # plt.plot(x, y, label="d=10mm")
    # plt.plot(smx[:500], smy[:500], label="d=1mm")
    # plt.plot(bigx[:500], bigy[:500], label="d=100mm")
    # plt.xlabel("x/d")
    # plt.ylabel("M")
    # plt.legend()
    
    # plt.plot(x,y-y)
    # plt.plot(smx[:500], 100*(smy[:500]-y[::2])/y[::2], label="d=1mm")
    # plt.plot(bigx[:500], 100*(bigy[:500]-y[::2])/y[::2], label="d=100mm")
    # plt.ylim(-1,1)
    # plt.xlim(-1,25)
    # plt.xlabel("x/d")
    # plt.ylabel("Mach error (%)")
    # plt.grid(True)
    
    
    # plt.show()
    # quit()
    
    # plt.plot(x, y)
    # plt.ylim(1,y[-1])
    # plt.xlim(0,x[-1])
    # plt.title("Centerline Mach number for PR=20")
    # plt.xlabel("x/D")
    # plt.ylabel("M")
    
    # plt.show()
    
    i_low = 10 # ~M=1.25
    M_low = 1.5
    i_low = np.argmin(np.abs(y-M_low))

    i_mid = 55 #~M=3
    M_mid = 2
    i_mid = np.argmin(np.abs(y-M_mid))
    print(x[i_mid])
    
    # i_high = i_mid

    M_high = 4
    i_high = np.argmin(np.abs(y-M_high))
    print(x[i_high])
    M_max=7
    i_max = np.argmin(np.abs(y-M_max))
    
    
    i_high
    # G = 1.40
    # f = lambda x, x0, A, C: A*(x-x0)**(G-1) - 0.5*(G+1)/(G-1)/(A*(x-x0)**(G-1))  +  C*(x-x0)**(-3*(G-1)) 

    
    # plt.plot(x, f(x, 0.40, 3.65, 0.20), label="Ashkenas & Sherman Model")
    # plt.ylim(0,max(y)*1.5)
    # plt.plot(x, y, label="Simulation")

    
    # f_var = farfield_model_var(x[i_high:], y[i_high:])
    # plt.plot(x, f_var(x), label="Modified A&S Model")
    # plt.legend()
    # plt.xlim(0,x[-1])
    # plt.xlabel("x/d")
    # plt.ylabel("M")      
    # plt.show()
    
    f_near = nearfield_model(x[:i_low], y[:i_low])
    # # plt.plot(x[:i_low],np.abs(f_near(x[:i_low])-y[:i_low])/y[:i_low])
    plt.plot(x[:i_low], f_near(x[:i_low]), 'r-', label="1+Ax+Bx^2+Cx^3")
    plt.plot(x[:i_low], y[:i_low], 'k--', lw=1)
    plt.xlabel("x/D")
    plt.ylabel("M")
    plt.legend()
    # plt.show()
    
    f_mixed = mixed_model(x[i_low:i_high],y[i_low:i_high])
    plt.plot(x[i_low:i_high],f_mixed(x[i_low:i_high]), 'g-', label="Ax+Bx^2+Cx^3+D")
    plt.plot(x[i_low:i_high], y[i_low:i_high], 'k--', lw=1)
    plt.xlabel("x/D")
    plt.ylabel("M")
    plt.legend()
    # plt.show()
        
    f_far = farfield_model_var(x[i_high:i_max], y[i_high:i_max])
    # plt.plot(x[i_mid:], np.abs(f_far(x[i_mid:])-y[i_mid:])/(y[i_mid:]))
    plt.plot(x[i_high:], f_far(x[i_high:]), label="Modified A&S Correlation")
    plt.plot(x[i_high:], y[i_high:], 'k--', lw=1, label = "Simulation")
    plt.legend()
    
    
    # plt.plot(x[i_low:i_high], y[i_low:i_high])
    
    
    
    plt.show()
    
    
    # f = lambda x, A, B, C, D: 1 + A*x**B + C*x**(3*B)
    # f2 = lambda x, x0, A, C, G: A*(x-x0)**(G-1) - 0.5*(G+1)/(G-1)/(A*(x-x0)**(G-1)) + C*(x-x0)**(-3*(G-1))
    
    # res = curve_fit(f, x, y, [1.5, 0.5, 1, 0])
    # # Model from the paper
    # res2 = curve_fit(f2, x, y, [-0.5, 0.5, 0.5, 1.4])
    
    # coeffs = res[0]
    # coeffs2 = res2[0]
    # print(coeffs2)
    # x_sam = np.linspace(0,4,len(y))
    # # plt.plot(x_sam, (f2(x_sam, *coeffs2)-y)/y)
    # # plt.plot([0,4], [0.01,0.01], 'r--')
    # # plt.plot([0,4], [-0.01,-0.01], 'r--')
    # plt.plot(x[:10], y[:10])
    # # plt.plot(xG, yG/yG[0])
    # plt.show()
    
    