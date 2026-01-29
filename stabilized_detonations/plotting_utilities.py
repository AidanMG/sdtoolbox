import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
from scipy.interpolate import interp1d


def extract_paraview_csv(filename, columns=None):
    data = pd.read_csv(filename)
    if columns:
        data = data[columns]
    return data

def plot_anim(df, x_var, y_var, time_col="Time", interval=50, repeat=True, x_lim=None, y_lim=None, title="", anim_filename=None):

    # ---- extract unique times ----
    times = np.sort(df[time_col].unique())

    # ---- pre-group for speed ----
    groups = dict(tuple(df.groupby(time_col)))

    # ---- global axis limits (prevents rescaling flicker) ----
    
    if x_lim:
        x_min, x_max = x_lim
    else:
        x_min = df[x_var].min()
        x_max = df[x_var].max() 
        
    if y_lim:
        y_min, y_max = y_lim
    else:
        y_min = df[y_var].min()
        y_max = df[y_var].max()        

    # ---- figure setup ----
    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=2)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.set_xlabel(x_var)
    ax.set_ylabel(y_var)
    ax_tit = ax.set_title(title)

    # ---- initialization ----
    def init():
        line.set_data([], [])
        if not title:
            ax_tit.set_text("")
        return line, ax_tit

    # ---- animation update ----
    def update(frame):
        t = times[frame]
        df_t = groups[t]

        x = df_t[x_var].to_numpy()
        y = df_t[y_var].to_numpy()

        # optional: sort for clean line plots
        order = np.argsort(x)
        x = x[order]
        y = y[order]

        line.set_data(x, y)
        if not title:
            ax_tit.set_text(f"{time_col} = {t:.6e}")

        return line, ax_tit

    anim = FuncAnimation(
        fig,
        update,
        frames=len(times),
        init_func=init,
        interval=interval,
        repeat=repeat,
        blit=True,
    )

    plt.show()
    if anim_filename:
        anim.save(anim_filename, writer="ffmpeg")

    return anim    

def average_last_timesteps(
    df,
    x_var,
    y_var,
    time_col="Time",
    n_last=10,
    n_grid=1000,
    kind="linear",
):
    """
    Time-average the last n_last timesteps while accounting for
    nonuniform and inconsistent spatial grids.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe.

    x_var : str
        Spatial coordinate (e.g. arc_length).

    y_var : str
        Variable to average (e.g. temperature, pressure).

    time_col : str
        Time column name.

    n_last : int
        Number of final timesteps to average.

    n_grid : int
        Number of points in reference interpolation grid.

    kind : str
        Interpolation type passed to scipy.interpolate.interp1d.

    Returns
    -------
    x_ref : ndarray
        Common arc-length grid.

    y_mean : ndarray
        Time-averaged variable on x_ref grid.

    y_std : ndarray
        Standard deviation (useful for unsteadiness).
    """

    # ---- sorted unique timesteps ----
    times = np.sort(df[time_col].unique())

    if n_last > len(times):
        raise ValueError("n_last exceeds available timesteps")

    selected_times = times[-n_last:]

    # ---- find global spatial bounds ----
    x_min = df[x_var].min()
    x_max = df[x_var].max()

    # ---- reference grid ----
    x_ref = np.linspace(x_min, x_max, n_grid)

    # ---- storage ----
    values = []

    # ---- loop over timesteps ----
    for t in selected_times:
        df_t = df[df[time_col] == t]

        x = df_t[x_var].to_numpy()
        y = df_t[y_var].to_numpy()

        # sort spatially
        order = np.argsort(x)
        x = x[order]
        y = y[order]

        # remove duplicates (important!)
        x_unique, idx = np.unique(x, return_index=True)
        y_unique = y[idx]

        # interpolate onto reference grid
        interp = interp1d(
            x_unique,
            y_unique,
            kind=kind,
            bounds_error=False,
            fill_value=np.nan,
        )

        y_interp = interp(x_ref)
        values.append(y_interp)

    values = np.array(values)

    # ---- average ignoring NaNs ----
    y_mean = np.nanmean(values, axis=0)
    y_std = np.nanstd(values, axis=0)

    return x_ref, y_mean, y_std


def get_M(data):
    data["M"] = np.sqrt(data["Velocity_u"]**2 + data["Velocity_v"]**2)/data["Speed_of_Sound"]
    return data

if __name__ == "__main__":
    data5 = extract_paraview_csv("PR5.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length"])
    data10 = extract_paraview_csv("PR10.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length"])
    data20 = extract_paraview_csv("PR20.csv", columns=["Time", "Velocity_u", "Velocity_v", "Speed_of_Sound", "arc_length"])
    data5 = get_M(data5)
    data10 = get_M(data10)
    data20 = get_M(data20)

    fig, ax = plt.subplots(2,2)

    data5_avg = average_last_timesteps(data5, "arc_length", "M", n_last=10)
    ax[0][0].set_xlim(0,2)
    ax[0][0].set_ylim(0,max(data5_avg[1]+data5_avg[2]))
    ax[0][0].set_title("PR5 averaged")
    ax[0][0].set_xlabel("arc_length")
    ax[0][0].set_ylabel('M')
    
    ax[0][0].plot(data5_avg[0], data5_avg[1])
    ax[0][0].plot(data5_avg[0], data5_avg[1]+data5_avg[2], 'b--', lw=0.5)
    ax[0][0].plot(data5_avg[0], data5_avg[1]-data5_avg[2], 'b--', lw=0.5)
    
    data10_avg = average_last_timesteps(data10, "arc_length", "M", n_last=10)
    ax[0][1].set_xlim(0,2)
    ax[0][1].set_ylim(0,max(data10_avg[1]+data10_avg[2]))
    ax[0][1].set_title("PR10 averaged")
    ax[0][1].set_xlabel("arc_length")
    ax[0][1].set_ylabel('M')
    
    ax[0][1].plot(data10_avg[0], data10_avg[1])
    ax[0][1].plot(data10_avg[0], data10_avg[1]+data10_avg[2], 'b--', lw=0.5)
    ax[0][1].plot(data10_avg[0], data10_avg[1]-data10_avg[2], 'b--', lw=0.5)

    data20_avg = average_last_timesteps(data20, "arc_length", "M", n_last=10, n_grid=2000)
    ax[1][0].set_xlim(0,4)
    ax[1][0].set_ylim(0,max(data20_avg[1]+data20_avg[2]))
    ax[1][0].set_title("PR20 averaged")
    ax[1][0].set_xlabel("arc_length")
    ax[1][0].set_ylabel('M')
    
    ax[1][0].plot(data20_avg[0], data20_avg[1])
    ax[1][0].plot(data20_avg[0], data20_avg[1]+data20_avg[2], 'b--', lw=0.5)
    ax[1][0].plot(data20_avg[0], data20_avg[1]-data20_avg[2], 'b--', lw=0.5)
    
    ax[1][1].set_xlabel("arc_length")
    ax[1][1].set_ylabel("M")
    ax[1][1].set_title("Comparative")
    ax[1][1].set_xlim(0,2)
    ax[1][1].set_ylim(0,max(data20_avg[1]))
    ax[1][1].plot(data5_avg[0], data5_avg[1])
    ax[1][1].plot(data10_avg[0], data10_avg[1])
    ax[1][1].plot(data20_avg[0], data20_avg[1])

    
    plt.show()
    
    # fig, ax = plt.subplots(2,2)
    # for a in ax.flatten():
    #     a.set_xlim(0,2)
    
    
    # ax[0][0].plot(data5_avg[0], data5_avg[1])
    # ax[0][0].plot(data20_avg[0], data20_avg[1], 'b--', lw=0.5)
    # ax[0][0].set_title("PR 5 vs. 20")
    # ax[0][0].set_ylabel("M")
    
    # ax[0][1].plot(data10_avg[0], data10_avg[1])
    # ax[0][1].plot(data20_avg[0], data20_avg[1], 'b--', lw=0.5)
    # ax[0][1].set_title("PR 10 vs. 20")
    
    # ax[1][0].plot(data5_avg[0], np.abs(data5_avg[1]-data20_avg[1][:1000])/np.abs(data20_avg[1][:1000]))
    # ax[1][1].plot(data10_avg[0], np.abs(data10_avg[1]-data20_avg[1][:1000])/np.abs(data20_avg[1][:1000]))
    # ax[1][0].set_yscale("log")
    # ax[1][1].set_yscale("log")    
    
    # ax[1][0].set_xlabel("arc length")
    # ax[1][0].set_ylabel("Relative Error")
    
    # ax[1][1].set_xlabel("arc length")
    
    # plt.show()
    
    # plot_anim(data5, "arc_length", "M", x_lim = (0, 0.75), interval=100, title="PR5", anim_filename="PR5.gif")
    # plot_anim(data10, "arc_length", "M", x_lim = (0, 1.5), interval=100, title="PR10", anim_filename="PR10.gif")
    # plot_anim(data20, "arc_length", "M", interval=100, title="PR20", anim_filename="PR20.gif")
    
    
