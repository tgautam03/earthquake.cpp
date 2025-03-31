# %%
import os
# Change the current working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from ipywidgets import interact, interactive, fixed, interact_manual
from matplotlib.colors import LinearSegmentedColormap
from IPython.display import display
import ipywidgets as widgets
import matplotlib.colors as mcolors

# Define a custom colormap (white to brown)
white_brown_cmap = LinearSegmentedColormap.from_list("white_brown", ["white", "brown"])

plt.style.use('dark_background')

# Enable interactive Matplotlib backend
%matplotlib widget

if __name__ == "__main__":
    # Load data
    x = np.frombuffer(open("../data/3d/x.bin", "rb").read(), dtype=np.float32)
    y = np.frombuffer(open("../data/3d/y.bin", "rb").read(), dtype=np.float32)
    z = np.frombuffer(open("../data/3d/z.bin", "rb").read(), dtype=np.float32)
    t = np.frombuffer(open("../data/3d/t.bin", "rb").read(), dtype=np.float32)

    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    dt = t[1] - t[0]

    # Grid Details
    nx = len(x)
    ny = len(y)
    nz = len(z)
    nt = len(t)

    c = np.frombuffer(open("../data/3d/c.bin", "rb").read(), dtype=np.float32)
    c = c.reshape(nz, ny, nx)

    # Solution
    sol = np.frombuffer(open("../data/3d/sol.bin", "rb").read(), dtype=np.float32)
    sol = sol.reshape(nt, nz, ny, nx)

    X, Z = np.meshgrid(x, z)
    zlims = max(np.abs(sol[:,:,0,:].max()), np.abs(sol[:,:,0,:].min()))

    # Create the initial plot
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    # Create a custom colormap
    cmap = mcolors.ListedColormap(['white', 'brown'])

    # Create a normalization that maps your data range to [0, 1]
    norm = mcolors.Normalize(vmin=c.min(), vmax=c.max())

    extent = [x[0], x[-1], y[-1], y[0]]

    im1 = axs[0].imshow(c[0], extent=extent, aspect='auto', cmap=cmap, norm=norm)
    src = axs[0].scatter((x[-1]-x[0])/2, (y[-1]-y[0])/1.25, color='red', marker='*', s=50, label="Explosive", zorder=2)
    rec1 = axs[0].scatter(x[0], y[0], color='green', marker='o', s=50, label="Location 1", zorder=2)
    rec2 = axs[0].scatter(x[0], y[0], color='purple', marker='o', s=50, label="Location 2", zorder=2)
    im2 = axs[0].imshow(sol[0,0,:,:], cmap='coolwarm', alpha=0.7, extent=extent, aspect='auto')

    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    title0 = axs[0].set_title("t: 0 sec")

    line1, = axs[1].plot([], [], linewidth=3, color='green')
    line2, = axs[1].plot([], [], linewidth=3, color='purple', linestyle="--")
    axs[1].set_xlim(0, t[-1])
    axs[1].set_ylabel("Pressure")
    axs[1].set_xlabel("Time (seconds)")
    title1 = axs[1].set_title("Pressure at a location")

    plt.tight_layout()

    def update_plot(x_val_g, y_val_g, x_val_p, y_val_p, z_val, t_val):
        x_i_g = int(x_val_g/dx)
        y_i_g = int(y_val_g/dy)
        x_i_p = int(x_val_p/dx)
        y_i_p = int(y_val_p/dy)
        z_i = int(z_val/dz)
        t_i = int(t_val/dt)

        # Update velocity field
        im1.set_data(c[z_i])

        # Update receiver location
        rec1.set_offsets([x_val_g, y_val_g])
        rec2.set_offsets([x_val_p, y_val_p])

        # Update wave field
        vlim = max(np.abs(sol[t_i,z_i,:,:].min()), np.abs(sol[t_i,z_i,:,:].max()))
        im2.set_data(sol[t_i,z_i,:,:])
        im2.set_clim(-vlim, vlim)

        # Update time series plot
        y_lims_g = max(np.abs(sol[:,z_i,y_i_g,x_i_g].max()), np.abs(sol[:,z_i,y_i_g,x_i_g].min()))
        y_lims_p = max(np.abs(sol[:,z_i,y_i_p,x_i_p].max()), np.abs(sol[:,z_i,y_i_p,x_i_p].min()))
        y_lims = max(y_lims_g, y_lims_p)
        line1.set_data(t[:t_i+1], sol[:t_i+1,z_i,y_i_g,x_i_g])
        line2.set_data(t[:t_i+1], sol[:t_i+1,z_i,y_i_p,x_i_p])
        axs[1].set_ylim(-1.1*y_lims, 1.1*y_lims)

        # Update titles
        title0.set_text(f"t: {t_val:.2f} sec")
        title1.set_text("Pressure at a location")

        fig.canvas.draw_idle()

    # Create an interactive slider and link it to the update function
    slider_x_g = widgets.SelectionSlider(options=x.tolist(), description="X Green (Km)")
    slider_y_g = widgets.SelectionSlider(options=y.tolist(), description="Y Green (Km)")
    slider_x_p = widgets.SelectionSlider(options=x.tolist(), description="X Purple (Km)")
    slider_y_p = widgets.SelectionSlider(options=y.tolist(), description="Y Purple (Km)")
    slider_z = widgets.SelectionSlider(options=z.tolist(), description="Z (Km)")
    slider_t = widgets.SelectionSlider(options=t[::10].tolist(), description="Time (s)")
    interactive_widget = widgets.interactive(update_plot, x_val_g=slider_x_g, y_val_g=slider_y_g, 
                                             x_val_p=slider_x_p, y_val_p=slider_y_p, 
                                             z_val=slider_z, t_val=slider_t)

    # Display the slider and initial plot
    display(interactive_widget)
    update_plot(x[0], y[0], x[0], y[0], z[0], t[0])  # Initialize with the first time steps