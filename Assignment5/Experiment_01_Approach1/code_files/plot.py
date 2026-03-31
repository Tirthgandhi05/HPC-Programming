# plot_exp01.py
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ===============================
# Configuration: grids and approaches
# ===============================
grids = [(250, 100), (500, 200), (1000, 400)]
approaches = ["immediate", "deferred"]
colors = {"immediate": "blue", "deferred": "red"}

# ===============================
# Loop over each grid and plot
# ===============================
for nx, ny in grids:
    plt.figure(figsize=(8,6))
    for approach in approaches:
        filename = f"exp01_{approach}_{nx}x{ny}.csv"
        try:
            df = pd.read_csv(filename)
            particles = df["Particles"].values
            total_time = df["TotalTime"].values

            plt.plot(particles, total_time, marker='o', color=colors[approach], label=approach.capitalize())
        except FileNotFoundError:
            print(f"Warning: {filename} not found. Skipping.")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Number of Particles (log scale)")
    plt.ylabel("Total Execution Time (s, log scale)")
    plt.title(f"Experiment 01: Execution Time vs Particle Count ({nx}x{ny} grid)")
    plt.legend()
    plt.grid(True, which="both", ls="--", lw=0.5)

    # Save figure
    plt.savefig(f"exp01_plot_{nx}x{ny}.png", dpi=300)
    print(f"Saved plot: exp01_plot_{nx}x{ny}.png")
    plt.close()
