import matplotlib.pyplot as plt
import pandas as pd
import os

grids = [(250, 100), (500, 200), (1000, 400)]
approaches = [("immediate", "Approach 2: Immediate Replacement"), 
              ("deferred", "Approach 1: Deferred Insertion")]

for nx, ny in grids:
    for approach_file, approach_title in approaches:
        filename = f"exp02_{approach_file}_{nx}x{ny}.csv"
        
        if not os.path.exists(filename):
            print(f"Skipping {filename}, not found.")
            continue
            
        df = pd.read_csv(filename)
        threads = df["Threads"]
        
        plt.figure(figsize=(7,5))
        
        # Plot Assignment 05 Speedup
        plt.plot(threads, df["Speedup"], marker='o', linewidth=2, color='blue', 
                 label=f"{approach_title} Speedup")
        
        # Plot Baseline (Assig 04) if it's the immediate file (to match assignment specs)
        if "BaselineSpeedup" in df.columns:
            plt.plot(threads, df["BaselineSpeedup"], marker='s', linewidth=2, color='orange', 
                     label="Baseline Mover (Assignment 04) Speedup")
        
        # Theoretical Max Line
        plt.plot(threads, threads, 'k--', label="Theoretical Maximum")
        
        plt.title(f"Mover Parallel Scalability ({nx}x{ny} Grid)")
        plt.xlabel("Number of OpenMP Threads")
        plt.ylabel("Speedup (S = T_serial / T_parallel)")
        plt.xticks(threads)
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        
        out_name = f"exp02_speedup_{approach_file}_{nx}x{ny}.png"
        plt.savefig(out_name, dpi=300)
        print(f"Generated {out_name}")
        plt.close()