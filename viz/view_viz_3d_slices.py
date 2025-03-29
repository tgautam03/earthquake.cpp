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

    # Solutions
    sol_vert = np.frombuffer(open("../data/3d/sol_vert.bin", "rb").read(), dtype=np.float32)
    sol_vert = sol_vert.reshape(nt, ny, nx)

    sol_surf = np.frombuffer(open("../data/3d/sol_surf.bin", "rb").read(), dtype=np.float32)
    sol_surf = sol_surf.reshape(nt, nz, nx)

    X, Z = np.meshgrid(x, z)
    X, Y = np.meshgrid(x, y)

    # Create the initial plot
    fig = plt.figure(figsize=(11, 5))
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.view_init(elev=25, azim=55)  # Elevation: 30°, Azimuth: 45°
    ax2 = fig.add_subplot(122)

    fig1, ax = plt.subplots(1, 1, figsize=(11, 4))

    # Initial surface plot
    def hide_3d_axes(ax):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax._axis3don = False  # This is the secret flag to hide 3D axes
    z_min = np.min(sol_surf)
    z_max = np.max(sol_surf)
    hide_3d_axes(ax1)
    surf = ax1.plot_surface(X, Z, sol_surf[0,:,:], cmap='coolwarm', alpha=0.7, vmin=z_min, vmax=z_max)
    rec1 = ax1.scatter(x[0], y[0], sol_surf[0,0,0], color='green', marker='o', s=50, label="Location 1", zorder=2)

    ax1.set_xlabel("x")
    ax1.set_ylabel("z")
    ax1.set_zlabel("u")
    ax1.set_xlabel("")
    ax1.set_ylabel("")
    ax1.set_zlabel("")
    ax1.set_zlim(z_min, z_max)
    title0 = ax1.set_title("t: 0 sec")

    ax2.axis('off')
    plt.tight_layout()

    # Initial vertical plot
    # Create a custom colormap
    cmap = mcolors.ListedColormap(['white', 'brown'])

    # Create a normalization that maps your data range to [0, 1]
    norm = mcolors.Normalize(vmin=c.min(), vmax=c.max())

    extent = [x[0], x[-1], y[-1], y[0]]

    im1 = ax2.imshow(c[0], extent=extent, aspect='auto', cmap=cmap, norm=norm)
    src = ax2.scatter((x[-1]-x[0])/2, (y[-1]-y[0])/1.25, color='red', marker='*', s=50, label="Earthquake", zorder=2)
    rec2 = ax2.scatter(x[0], y[0], color='purple', marker='o', s=50, label="Location 2", zorder=2)
    rec0 = ax2.scatter(x[0], 0, color='green', marker='o', s=50, label="Location 1", zorder=2)
    im2 = ax2.imshow(sol_vert[0,:,:], cmap='coolwarm', alpha=0.7, extent=extent, aspect='auto')

    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    title1 = ax2.set_title(f"t: {0} sec")

    line1, = ax.plot([], [], linewidth=3, color='green')
    line2, = ax.plot([], [], linewidth=3, color='purple', linestyle="--")
    ax.set_xlim(0, t[-1])
    ax.set_ylabel("u")
    ax.set_xlabel("Time (seconds)")
    title2 = ax.set_title(f"t: {0} sec")

    plt.tight_layout()

    def update_plot(x_val_g, z_val_g, x_val_p, y_val_p, t_val):
        global rec1  # Declare rec1 as global

        x_i_g = int(x_val_g/dx)
        z_i_g = int(z_val_g/dx)
        x_i_p = int(x_val_p/dx)
        y_i_p = int(y_val_p/dy)
        t_i = int(t_val/dt)

        # Update receiver location
        rec1._offsets3d = ([x_val_g], [0], [sol_surf[t_i,z_i_g,x_i_g]])

        # Update surface plot
        ax1.clear()
        hide_3d_axes(ax1)  # Re-apply axis hiding after clear()
        surf = ax1.plot_surface(X, Z, sol_surf[t_i,:,:], cmap='coolwarm', alpha=0.7, vmin=z_min, vmax=z_max)
        rec1 = ax1.scatter(x_val_g, z_val_g, sol_surf[t_i,z_i_g,x_i_g], color='green', marker='o', s=50, label="Location 1", zorder=2)
        ax1.set_xlabel("x")
        ax1.set_ylabel("z")
        ax1.set_zlabel("u")
        ax1.set_zlim(z_min, z_max)

        # Update receiver location
        rec0.set_offsets([x_val_g, 0])
        rec2.set_offsets([x_val_p, y_val_p])

        # Update wave field
        vlim = max(np.abs(sol_vert[t_i,:,:].min()), np.abs(sol_vert[t_i,:,:].max()))
        im2.set_data(sol_vert[t_i,:,:])
        im2.set_clim(-vlim, vlim)

        # Update time series plot
        y_lims_g = max(np.abs(sol_surf[:,z_i_g,x_i_g].max()), np.abs(sol_surf[:,z_i_g,x_i_g].min()))
        y_lims_p = max(np.abs(sol_vert[:,y_i_p,x_i_p].max()), np.abs(sol_vert[:,y_i_p,x_i_p].min()))
        y_lims = max(y_lims_g, y_lims_p)
        line1.set_data(t[:t_i+1], sol_surf[:t_i+1,z_i_g,x_i_g])
        line2.set_data(t[:t_i+1], sol_vert[:t_i+1,y_i_p,x_i_p])
        ax.set_ylim(-1.1*y_lims, 1.1*y_lims)

        # Update titles
        title0.set_text(f"t: {t_val:.2f} sec")
        title1.set_text(f"t: {t_val:.2f} sec")
        title2.set_text(f"t: {t_val:.2f} sec")

        fig.canvas.draw_idle()

    # Create an interactive slider and link it to the update function
    slider_x_g = widgets.SelectionSlider(options=x.tolist(), description="X Green (Km)")
    slider_z_g = widgets.SelectionSlider(options=z.tolist(), description="Z Green (Km)")
    slider_x_p = widgets.SelectionSlider(options=x.tolist(), description="X Purple (Km)")
    slider_y_p = widgets.SelectionSlider(options=y.tolist(), description="Y Purple (Km)")
    slider_t = widgets.SelectionSlider(options=t.tolist(), description="Time (s)")
    interactive_widget = widgets.interactive(update_plot, x_val_g=slider_x_g, z_val_g=slider_z_g, x_val_p=slider_x_p, y_val_p=slider_y_p, t_val=slider_t)

    # Display the slider and initial plot
    display(interactive_widget)
    update_plot(x[0], z[0], x[0], y[0], t[0])  # Initialize with the first time steps