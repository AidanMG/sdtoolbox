import pandas as pd
import matplotlib.pyplot as plt
import cantera as ct

# --- Load CSVs ---
csv1 = "PR2_data.csv"
csv2 = "PR5_data.csv"

df1 = pd.read_csv(csv1)
df2 = pd.read_csv(csv2)

# --- Extract ---
x1 = df1["Points:0"]
x2 = df2["Points:0"]

p1, p2   = df1["Pressure"],    df2["Pressure"]
rho1,rho2 = df1["Density"],    df2["Density"]
T1, T2   = df1["Temperature"], df2["Temperature"]
u1, u2   = df1["Velocity_u"],    df2["Velocity_u"]

# --- Create 2×2 subplot grid ---
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

gas = ct.Solution('cti/Lietal_2003.yaml')

def Mach(u, T):
    R = ct.gas_constant / gas.mean_molecular_weight
    gamma = gas.cp / gas.cv
    a = (0.1*gamma * R * T)**0.5
    return u / a


# Pressure
axes[0,0].plot(x1, p1/p1[0], label=csv1)
axes[0,0].plot(x2, p2/p2[0], label=csv2)
axes[0,0].set_title("P/Pe vs Location")
axes[0,0].set_xlabel("X")
axes[0,0].set_ylabel("Pressure")
axes[0,0].grid(True)
axes[0,0].legend()

# Density
axes[0,1].plot(x1, rho1/rho1[0], label=csv1)
axes[0,1].plot(x2, rho2/rho2[0], label=csv2)
axes[0,1].set_title("rho/rho_e vs Location")
axes[0,1].set_xlabel("X")
axes[0,1].set_ylabel("Density")
axes[0,1].grid(True)
axes[0,1].legend()

# Temperature
axes[1,0].plot(x1, T1/T1[0], label=csv1)
axes[1,0].plot(x2, T2/T2[0], label=csv2)
axes[1,0].set_title("T/T_e vs Location")
axes[1,0].set_xlabel("X")
axes[1,0].set_ylabel("Temperature")
axes[1,0].grid(True)
axes[1,0].legend()

# Velocity
axes[1,1].plot(x1, u1/u1[0], label=csv1)
axes[1,1].plot(x2, u2/u2[0], label=csv2)
axes[1,1].set_title("X-velocity vs Location")
axes[1,1].set_xlabel("X")
axes[1,1].set_ylabel("Velocity")
axes[1,1].grid(True)
axes[1,1].legend()

plt.tight_layout()
plt.show()
