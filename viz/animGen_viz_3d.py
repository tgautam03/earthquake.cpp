import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from moviepy import *
from tqdm import tqdm
import math

# Change the current working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

pio.kaleido.scope.default_format = "png"
pio.kaleido.scope.default_width = 480
pio.kaleido.scope.default_height = 480
pio.kaleido.scope.default_scale = 0.8

def create_frame(t_i, num_iterations):
    angle = 2 * math.pi * t_i / num_iterations
    camera_x = 2 * math.cos(angle)
    camera_y = 0.5
    camera_z = 2 * math.sin(angle)
    
    fig = go.Figure(data=go.Volume(
        x=X_flat,
        y=Y_flat,
        z=Z_flat,
        value=sol_needed[t_i],
        isomin=sol_needed[t_i].min(),
        isomax=sol_needed[t_i].max(),
        opacity=0.1,
        surface_count=10,
    ))
    fig.update_layout(
        scene_camera=dict(
            eye=dict(x=camera_x, y=camera_y, z=camera_z),
            center=dict(x=0, y=0, z=0),
            up=dict(x=0, y=0, z=1)
        ),
        scene=dict(
            aspectmode='cube'
        )
    )
    fig.update_layout(constant_layout)
    pio.write_image(fig, f"./frames_3d/frame_{t_i}.png")
    fig.data = []
    fig = None

def create_frame_wrapper(args):
    t_i, num_iterations = args
    create_frame(t_i, num_iterations)

if __name__ == "__main__":
    x = np.frombuffer(open("../data/3d/x.bin", "rb").read(), dtype=np.float32)
    y = np.frombuffer(open("../data/3d/y.bin", "rb").read(), dtype=np.float32)
    z = np.frombuffer(open("../data/3d/z.bin", "rb").read(), dtype=np.float32)
    t = np.frombuffer(open("../data/3d/t.bin", "rb").read(), dtype=np.float32)

    nx, ny, nz, nt = len(x), len(y), len(z), len(t)

    c = np.frombuffer(open("../data/3d/c.bin", "rb").read(), dtype=np.float32).reshape(nz,ny,nx)

    sol = np.frombuffer(open("../data/3d/sol.bin", "rb").read(), dtype=np.float32).reshape(nt, nz, ny, nx)
    step = 10
    sol_needed = sol[::step]
    del sol

    print("===============")
    print("=== Details ===")
    print("===============")
    print(f"x:({x[0].item(), x[-1].item()}), y:({y[0].item(), y[-1].item()}), z:({z[0].item(), z[-1].item()}), t:({t[0].item(), t[-1].item()})")
    print(f"nx:{sol_needed.shape[-1]}, ny:{sol_needed.shape[-2]}, nz:{sol_needed.shape[-3]}, nt:{nt}")
    print(f"Solution shape: {sol_needed.shape}")

    X, Y, Z = np.mgrid[:sol_needed.shape[-1], :sol_needed.shape[-2], :sol_needed.shape[-3]]

    X_flat, Y_flat, Z_flat = X.flatten(), Y.flatten(), Z.flatten()

    constant_layout = dict(
        paper_bgcolor="black", plot_bgcolor="black",
        scene=dict(
            xaxis=dict(showticklabels=False, showaxeslabels=False, showline=False, zeroline=True, backgroundcolor="rgb(0, 0, 0)", gridcolor="black", showbackground=True, zerolinecolor="black"),
            yaxis=dict(showticklabels=False, showaxeslabels=False, showline=False, zeroline=True, backgroundcolor="rgb(0, 0, 0)", gridcolor="black", showbackground=True, zerolinecolor="black"),
            zaxis=dict(showticklabels=False, showaxeslabels=False, showline=False, zeroline=True, backgroundcolor="rgb(0, 0, 0)", gridcolor="black", showbackground=True, zerolinecolor="black"),
        )
    )
    
    num_iterations = len(sol_needed)
    sol_needed = sol_needed.reshape(num_iterations, -1)

    print("Generating frames_3d (this may take a while)...")

    for t_i in tqdm(range(216, num_iterations)):
        create_frame(t_i, num_iterations)

    print(f"All {num_iterations} Frames generated!")

    print("Creating mp4 from frames_3d...")
    image_folder = './frames_3d/'
    fps = 30

    def sorted_alphanumeric(data):
        import re
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(data, key=alphanum_key)

    image_files = [os.path.join(image_folder, img) 
                   for img in sorted_alphanumeric(os.listdir(image_folder))
                   if img.endswith(".png")]

    clip = ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile("wave_3d.mp4")
    print("mp4 generated: wave_3d.mp4")