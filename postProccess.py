import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Directory containing the trajectory files
trajectory_dir = "build/Trajectories/ExerciseC/3/ST_10/data/"

# Domain boundaries
x_min, x_max = 0, 0.00385802
y_min, y_max = 0, 0.000771605

def load_trajectory(file_path):
    """Load a single trajectory file."""
    return pd.read_csv(file_path, delim_whitespace=True, comment="#", header=None,
                       names=["t", "x", "y", "z", "u", "v", "m", "T"])

# Initialize empty numpy arrays for storing the coordinates and temperatures
x_coords_all = []
y_coords_all = []
temps_all = []

# Load all trajectory files in the directory
for file_name in os.listdir(trajectory_dir):
    file_path = os.path.join(trajectory_dir, file_name)
    if os.path.isfile(file_path):  # Ensure it's a file
        data = load_trajectory(file_path)
        x_coords_all.append(data["x"].values)
        y_coords_all.append(data["y"].values)
        temps_all.append(data["T"].values)

# Convert lists to single arrays for plotting and KDE
x_coords_all = np.concatenate(x_coords_all)
y_coords_all = np.concatenate(y_coords_all)
temps_all = np.concatenate(temps_all)

# Gaussian KDE for temperature field
points = np.vstack([x_coords_all, y_coords_all])
values = temps_all
kde = gaussian_kde(points, weights=values)

# Create a grid for the temperature field
x_grid = np.linspace(x_min, x_max, 200)
y_grid = np.linspace(y_min, y_max, 200)
x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)
positions = np.vstack([x_mesh.ravel(), y_mesh.ravel()])

# Evaluate the KDE
temperature_field = kde(positions).reshape(x_mesh.shape)

# Normalize the temperature field and scale
temperature_field_normalized = (temperature_field / np.max(temperature_field)) * 370

# Plot the scalar field and scatter points
plt.figure(figsize=(12, 8))
plt.title("Normalized Temperature Field with Gaussian KDE and Scatter Plot")
plt.xlabel("x")
plt.ylabel("y")

# Plot the normalized temperature field
plt.pcolormesh(x_mesh, y_mesh, temperature_field_normalized, shading="auto", cmap="cool")
plt.colorbar(label="Normalized Temperature Field (K)")

# Plot the original scatter points
scatter = plt.scatter(x_coords_all, y_coords_all, c=temps_all, cmap="cool", edgecolor="k", s=10, label="Original Data")
plt.colorbar(scatter, label="Original Temperature (K)")
plt.xlim((x_min, x_max))
plt.ylim((y_min, y_max))
plt.legend()
plt.grid()
plt.show()
