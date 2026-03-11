import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load CSV
df = pd.read_csv("exp1_results.csv")

configs = df["Config"].unique()

for config in configs:

    subset = df[df["Config"] == config]

    particles = subset["Particles"]
    time = subset["TotalInterpolationTime"]

    plt.figure()
    plt.loglog(particles, time, marker='o')

    plt.xlabel("Number of Particles (log scale)")
    plt.ylabel("Total Interpolation Time (log scale)")
    plt.title(f"Experiment 01 - Configuration {config}")
    plt.grid(True)

    plt.savefig(f"exp1_config_{config}.png")
    plt.close()

print("Plots saved.")
