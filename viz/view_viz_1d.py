# %%
import os
# Change the current working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import display
import matplotlib.colors as mcolors

plt.style.use('dark_background')

# Enable interactive Matplotlib backend
%matplotlib widget

MY_GOLD = "#DAA520"

if __name__ == "__main__":
    FILE_EXT = "3pt_00"
    # FILE_EXT = "3pt_da"
    # FILE_EXT = "5pt_da"
    # FILE_EXT = "7pt_da"
    
    y = np.frombuffer(open("../data/1d/y_{}.bin".format(FILE_EXT), "rb").read(), dtype=np.float32)
    dy = y[1] - y[0]
    t = np.frombuffer(open("../data/1d/t_{}.bin".format(FILE_EXT), "rb").read(), dtype=np.float32)
    dt = t[1] - t[0]


    c = np.frombuffer(open("../data/1d/c_{}.bin".format(FILE_EXT), "rb").read(), dtype=np.float32)

    # Grid Details
    nx = len(y)
    nt = len(t)

    # Solution
    sol = np.frombuffer(open("../data/1d/sol_{}.bin".format(FILE_EXT), "rb").read(), dtype=np.float32)
    sol = sol.reshape(nt, nx)
    
    y_lims = max(np.abs(sol.max()), np.abs(sol.min()))

    # Create the figure and axes once
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    # Create a custom colormap
    cmap = mcolors.ListedColormap(['white', 'brown'])

    # Create a normalization that maps your data range to [0, 1]
    norm = mcolors.Normalize(vmin=c.min(), vmax=c.max())

    # Pre-create the plots
    line1, = axs[0].plot([], [], linewidth=5, color=MY_GOLD, zorder=1)  # Line plot for u(y,t) vs spatial grid points
    im = axs[0].imshow(np.ones((nx,nx))*np.expand_dims(c, axis=-1), extent=[-1.1*y_lims, 1.1*y_lims, y.max(), y.min()], 
                       aspect='auto', cmap=cmap, norm=norm)
    src = axs[0].scatter(0, int((y[-1]-y[0])/2), color='red', marker='*', s=50, label="Explosive", zorder=2)
    rec1 = axs[0].scatter(0, 0, color='green', marker='o', s=50, label="Location 1", zorder=2)
    rec2 = axs[0].scatter(0, 0, color='purple', marker='o', s=50, label="Location 2", zorder=2)
    axs[0].set_xlim(-1.1 * y_lims, 1.1 * y_lims)
    axs[0].set_ylim(y.min(), y.max())  # Set limits based on spatial grid points
    axs[0].set_ylabel("y")
    axs[0].set_xlabel("u(y,t)")
    axs[0].invert_yaxis()
    # cbar = fig.colorbar(im, ax=axs)
    # cbar.set_label("Velocity (km/s)")

    line2, = axs[1].plot([], [], color='green', linewidth=3)  # Line plot for u(y,t) at location 1 over time
    line3, = axs[1].plot([], [], color='purple', linewidth=3, linestyle="--")  # Line plot for u(y,t) at location 2 over time
    axs[1].set_xlim(0, t[-1])
    axs[1].set_ylim(-1.1 * y_lims, 1.1 * y_lims)
    axs[1].set_ylabel("u(y,t)")
    axs[1].set_xlabel("Time (seconds)")

    plt.tight_layout()

    # Function to update the plots dynamically
    def update_plot(t_val, loc1, loc2):
        t_i = int(t_val / dt)
        loc1_i = int(loc1 / dy)
        loc2_i = int(loc2 / dy)
        
        # Update line plot data for u(y,t) vs spatial grid points
        line1.set_data(sol[t_i], y)  # Update u(y,t) (sol[t_i]) vs spatial grid points (y)
        
        # Update line plot data for u(y,t) under the surface over time
        line2.set_data(t[:t_i + 1], sol[:t_i + 1, loc1_i])  # Update time series plot

        # Update line plot data for u(y,t) over surface over time
        line3.set_data(t[:t_i + 1], sol[:t_i + 1, loc2_i])  # Update time series plot

        # Update positions of rec1 and rec2
        rec1.set_offsets([0, y[loc1_i]])
        rec2.set_offsets([0, y[loc2_i]])
        
        # Update titles dynamically
        axs[0].set_title(f"u(y,t) Distribution at t: {t_val:.2f} sec")
        # axs[1].set_title(f"u(y,t) at location 1 over Time (t={t_val:.2f})")
        # axs[2].set_title(f"u(y,t) at location 2 over Time (t={t_val:.2f})")
        
        fig.canvas.draw_idle()  # Redraw only updated elements

    # Create an interactive slider and link it to the update function
    slider_t = widgets.SelectionSlider(options=t.tolist(), description="Time (s)")
    slider_loc1 = widgets.SelectionSlider(options=y.tolist(), description="Loc 1 (Km)")
    slider_loc2 = widgets.SelectionSlider(options=y.tolist(), description="Loc 2 (Km)")
    interactive_widget = widgets.interactive(update_plot, t_val=slider_t, loc1=slider_loc1, loc2=slider_loc2)

    # Display the slider and initial plot
    display(interactive_widget)
    update_plot(t[0], 0, 0)  # Initialize with the first time step
# %%
