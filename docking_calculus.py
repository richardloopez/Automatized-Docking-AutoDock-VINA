#!/usr/bin/env python3

# Author: Cristina Garcia Iriepa
# Citation: If you use this code, please cite Garcia-Iriepa, C

import os
import subprocess
import time

def run_vina():
    """
    Run AutoDock Vina 50 times and save the output in a specific folder.
    """
    folder = "receptor_ligand"
    
    # Create the folder
    os.makedirs(folder, exist_ok=True)

    for i in range(1, 51):  # Run 50 times, from 1 to 50
        output_file = f"./{folder}/{folder}_{i}.pdbqt"
        log_file = f"./{folder}/{folder}_{i}.txt"

        # Run Vina command
        command = f"vina --config config.txt --out {output_file}"

        try:
            # Run the command and capture the output
            with open(log_file, 'w') as log:
                subprocess.run(command, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
            print(f"Completed run {i}/50")
        except subprocess.CalledProcessError as e:
            print(f"Error in run {i}: {e}")

        # Wait for 10 seconds before the next iteration
        time.sleep(10)

if __name__ == "__main__":
    run_vina()

print("All Vina runs completed.")
