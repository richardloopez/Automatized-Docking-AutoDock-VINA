#!/usr/bin/env python3

# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
# Citation: If you use this code, please cite Lopez-Corbalan, R

import os
import shutil
import re

def extract_energy(file_path):
    """
    Extract the energy value from the second line of a PDBQT file.

    Args:
        file_path (str): Path to the PDBQT file.

    Returns:
        str: Extracted energy value without decimal point.
    """
    with open(file_path, 'r') as file:
        # Skip the first line
        next(file)
        # Read the second line
        second_line = next(file)
        # Extract the energy value
        energy_match = re.search(r'-?\d+\.\d+', second_line)
        if energy_match:
            # Remove the decimal point and return
            return energy_match.group().replace('.', '')
    return None

def sort_poses():
    """
    Sort PDBQT files based on their energy values and move them to a new directory.
    """
    # Name of the output directory
    output_dir = "energy_sorted"

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all .pdbqt files in the current directory
    for pdbqt_file in os.listdir('.'):
        if pdbqt_file.endswith('.pdbqt'):
            # Extract energy from the file
            energy = extract_energy(pdbqt_file)
            if energy is not None:
                # Create a new file name with the energy
                new_file_name = f"E{energy}_P{pdbqt_file}"

                # Move and rename the file to the new directory
                shutil.move(pdbqt_file, os.path.join(output_dir, new_file_name))
                print(f"Moved and renamed: {pdbqt_file} -> {new_file_name}")
            else:
                print(f"Could not extract energy from {pdbqt_file}")

if __name__ == "__main__":
    sort_poses()

