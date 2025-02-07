#!/usr/bin/env python3

# Author: Lorenzo Gramolini
# GitHub: github.com/lorenzo-gram
# Citation: If you use this code, please cite Gramolini, L

import os
import re

def add_filename_to_model(input_file, output_file):
    """
    Add the filename to each MODEL line in a PDBQT file.

    Args:
        input_file (str): Path to the input PDBQT file.
        output_file (str): Path to the output PDBQT file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('MODEL'):
                # Add filename to the MODEL line
                outfile.write(f"{line.strip()} {os.path.basename(input_file)}\n")
            else:
                # Write other lines as they are
                outfile.write(line)

def process_pdbqt_files():
    """
    Process all PDBQT files in the current directory,
    adding the filename to each MODEL line.
    """
    for filename in os.listdir('.'):
        if filename.endswith('.pdbqt'):
            temp_file = f"temp_{filename}"
            add_filename_to_model(filename, temp_file)
            os.replace(temp_file, filename)
            print(f"Processed: {filename}")

if __name__ == "__main__":
    process_pdbqt_files()

