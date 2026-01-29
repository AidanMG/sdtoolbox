import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -----------------------------
# USER SETTINGS
# -----------------------------
csv_file = "centerline.csv"

arc_col = "arc_length"

# Four quantities to plot (must match column names in CSV)
quantities = [
    "Pressure",
    "Temperature",
    "M",
    
    # "Y_03__H2O"
]

# Identity is the default if not specified
def identity(x):
    return x

transforms = {
    "Pressure": lambda y: y / 101325.0,   # normalize by 1 atm
    "Temperature":  lambda y: y / 300,      # normalize by 300K
    "M": identity,
    # "velocity": lambda y: y / 300.0
}

# Which subplots should be log-scaled
log_scale = {
    "Pressure": False,
    "Temperature": False,
    "M": False,
    # "Y_03__H2O": True
}

# Animation parameters
interval_ms = 500        # time between frames
repeat = True

# -----------------------------
# LOAD CSV
# -----------------------------
print("LOADING CSV...")
data = np.genfromtxt(csv_file, delimiter=",", names=True)
print(data.dtype.names)
arc = data[arc_col]

# Detect timestep boundaries (arc_length ~ 0)
# Use tolerance to avoid floating-point issues
tol = 1e-12
timestep_starts = np.where(np.abs(arc) < tol)[0]

# Add end index for slicing
timestep_starts = np.append(timestep_starts, len(arc))

n_timesteps = len(timestep_starts) - 1

# -----------------------------
# PRE-SPLIT DATA BY TIMESTEP
# -----------------------------
print("SPLITTING DATA...")
timesteps = []

for i in range(n_timesteps):
    i0 = timestep_starts[i]
    i1 = timestep_starts[i + 1]

    step = {
        "arc_length": arc[i0:i1]
    }
    
    for q in quantities:
        step[q] = data[q][i0:i1]
    timesteps.append(step)

# -----------------------------
# SET UP FIGURE
# -----------------------------
print("PLOTTING...")
y_limits = {}

for q in quantities:
    y = data[q]
    
    # Apply transform
    y = transforms.get(q, identity)(y)    

    # Optional: guard against NaNs / infs
    y = y[np.isfinite(y)]

    ymin = y.min()
    ymax = y.max()

    # Optional padding (5%)
    pad = 0.05 * (ymax - ymin) if ymax > ymin else 1.0
    y_limits[q] = (ymin - pad, ymax + pad)

fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
axes = axes.flatten()

lines = []

for ax, q in zip(axes, quantities):
    line, = ax.plot([], [], lw=2)
    lines.append(line)

    ax.set_ylabel(q)
    ax.grid(True)

    if log_scale.get(q, False):
        ax.set_yscale("log")

        # Ensure limits are valid for log scale
        ymin, ymax = y_limits[q]
        print(ymin,ymax)
        ax.set_ylim(ymin, ymax)
    else:
        ax.set_ylim(*y_limits[q])

axes[2].set_xlabel("arc_length")
axes[3].set_xlabel("arc_length")

title = fig.suptitle("")

# -----------------------------
# INITIALIZATION FUNCTION
# -----------------------------
def init():
    for line in lines:
        line.set_data([], [])
    title.set_text("")
    return lines + [title]

# -----------------------------
# ANIMATION UPDATE FUNCTION
# -----------------------------
def update(frame):
    step = timesteps[frame]

    x = step["arc_length"]

    for line, q in zip(lines, quantities):
        y = step[q]
        y=transforms.get(q,identity)(y)
        
        y = np.clip(y, 1e-12, None)
        line.set_data(x, y)

    # Update axis limits dynamically
    for ax, q in zip(axes, quantities):
        ax.relim()
        ax.autoscale_view()

    title.set_text(f"Timestep {frame + 1} / {n_timesteps}")

    return lines + [title]

# -----------------------------
# CREATE ANIMATION
# -----------------------------
print("WRITING ANIMATION...")
ani = FuncAnimation(
    fig,
    update,
    frames=n_timesteps,
    init_func=init,
    interval=interval_ms,
    blit=False,
    repeat=repeat
)
from matplotlib.animation import PillowWriter

ani.save(
    "animation.gif",
    writer=PillowWriter(fps=1000/interval_ms),
    dpi=150
)


plt.tight_layout()
plt.show()