# Automated Docking Workflow with AutoDock Vina

This repository contains a set of Python scripts to automate the docking process using AutoDock Vina. The workflow includes preparing input files, running multiple docking simulations, analyzing the results based on RMSD, and selecting the most promising poses based on energy.

## Table of Contents

1.  [Overview](#overview)
2.  [Requirements](#requirements)
3.  [Installation](#installation)
4.  [Usage](#usage)
5.  [File Descriptions](#file-descriptions)
6.  [Configuration](#configuration)
7.  [Workflow](#workflow)
8.  [Authors](#authors)
9.  [Citation](#citation)

## Overview

This workflow streamlines the process of running multiple AutoDock Vina simulations, filtering poses based on RMSD, and selecting the best poses based on their energy scores.  The scripts handle tasks such as:

*   Running AutoDock Vina multiple times.
*   Adding the filename to each MODEL line in the output `.pdbqt` files.
*   Calculating RMSD values between poses to identify unique binding modes.
*   Merging selected poses into a single file.
*   Sorting poses by energy.

## Requirements

To use this workflow, you need the following:

*   Python 3.x
*   AutoDock Vina
*   Python libraries: `os`, `subprocess`, `glob`, `shutil`, `re`, `time` (most are standard libraries)
*   `numpy` (if you intend to add any numpy-related functionality)

## Installation

1.  Clone the repository:

    ```
    git clone <repository_url>
    cd <repository_name>
    ```

2.  Ensure AutoDock Vina is installed and accessible in your system's PATH.
3.  Install any missing Python libraries (though most are standard):

    ```
    pip install numpy  # If needed
    ```

## Usage

1.  Prepare your ligand and receptor files in `.pdbqt` format.
2.  Create a directory for your docking experiment.
3.  Place the ligand `.pdbqt` file, the receptor `.pdbqt` file, all the Python scripts (`Automated_Docking_VINA.py`, `docking_calculus.py`, `add_name.py`, `analysis_docking_general_auto.py`, `analysis_docking_particular_auto.py`, `sort_poses.py`), and the `config.txt` file into this directory.
4.  Modify the `config.txt` file and the scripts as needed to match your specific docking parameters (see [Configuration](#configuration)).
5.  Run the `Automated_Docking_VINA.py` script using `sbatch` or a similar job scheduler:

    ```
    sbatch Automated_Docking_VINA.py
    ```

## File Descriptions

*   `Automated_Docking_VINA.py`:  The main script that controls the entire docking workflow. It launches `docking_calculus.py`, copies scripts, executes pose analysis, and sorts poses.
*   `docking_calculus.py`: Runs AutoDock Vina multiple times (default: 50) based on parameters in `config.txt`. Creates a folder (default: `receptor_ligand`) to store output files.
*   `add_name.py`: Adds the filename to each `MODEL` line in the `.pdbqt` files. This is crucial for tracking the origin of each pose after merging.
*   `analysis_docking_general_auto.py`: Filters poses based on RMSD values, keeping only diverse poses, and merges the selected poses into a `merged.pdbqt` file.
*   `analysis_docking_particular_auto.py`: Performs further analysis on the merged poses to prepare files for final selection.
*   `sort_poses.py`: Sorts the final poses based on their energy scores and places the selected poses into a final directory.
*   `config.txt`: Configuration file for AutoDock Vina, specifying the receptor, ligand, search box coordinates and dimensions, and number of modes.

## Configuration

### `config.txt`

This file contains the core docking parameters:

*   `receptor`: Path to the receptor `.pdbqt` file.
*   `ligand`: Path to the ligand `.pdbqt` file.
*   `center_x`, `center_y`, `center_z`: Coordinates for the center of the search box.
*   `size_x`, `size_y`, `size_z`: Dimensions of the search box.
*   `num_modes`: Number of poses to generate per docking run.

### Adaptable Parameters in Scripts

| File                                  | Parameter                               | Description                                                                                                                          |
| ------------------------------------- | --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `config.txt`                          | Ligand and receptor names               | Specify the names (and paths, if necessary) of the ligand and receptor `.pdbqt` files.                                            |
| `config.txt`                          | Box parameters                          | Define the coordinates and dimensions of the search box.  Crucial for targeting the correct binding site.                             |
| `config.txt`                          | `num_modes`                             | Set the number of poses generated by Vina for each individual docking run.                                                          |
| `docking_calculus.py`                 | `folder`                                | Name of the folder where the individual docking run results will be stored.                                                        |
| `docking_calculus.py`                 | Number of docking runs                  | Modify the `range` in the `for` loop to change the number of times Vina is executed.                                                |
| `analysis_docking_general_auto.py`    | `pose_inicial`, `pose_final`            | Define the range of poses (ranked by energy) that will be considered for RMSD analysis and merging.                                  |
| `analysis_docking_general_auto.py`    | RMSD threshold                          | Adjust the RMSD value (currently set to `< 2`) used to determine if two poses are considered similar.  Lower values = more diversity. |

## Workflow

1.  **Preparation**:
    *   Create a dedicated folder for your docking project.
    *   Place all `.py` scripts, the `config.txt` file, and the ligand and receptor `.pdbqt` files inside.
2.  **Docking Calculations**:
    *   Run `Automated_Docking_VINA.py`. This script:
        *   Launches `docking_calculus.py`, which performs multiple Vina docking runs.
        *   Creates a subfolder (e.g., `receptor_ligand`) to store the results of each run.
3.  **Analysis and Merging**:
    *   `Automated_Docking_VINA.py` then:
        *   Copies the analysis scripts into the results subfolder.
        *   Executes `add_name.py` to add filenames to the `MODEL` lines of the `.pdbqt` files.
        *   Runs `analysis_docking_general_auto.py` to filter poses based on energy and RMSD, merging the selected diverse poses into `merged.pdbqt`.
4.  **Refinement and Selection**:
    *   The script creates a `final_selection_1` folder and moves `analysis_docking_particular_auto.py`, `sort_poses.py`, and `merged.pdbqt` into it.
    *   `analysis_docking_particular_auto.py` is executed for further processing.
5.  **Sorting and Output**:
    *   An `ultimate_selection` folder is created. Relevant `.pdbqt` files are moved into this folder.
    *   `sort_poses.py` is executed to sort the poses based on their energy scores. The final, sorted results are placed in the `resultado_energia` folder.
6.  **Visualization and Analysis**:
    *   Download the files from the `ultimate_selection` and/or `resultado_energia` folder.
    *   Analyze the poses using Discovery Studio, PyMOL, or your preferred molecular visualization software.

## Authors

*   Richard Lopez Corbalan (github.com/richardloopez)
*   Cristina Garcia Iriepa
*   Lorenzo Gramolini (github.com/lorenzo-gram)

## Citation

If you use this workflow in your research, please cite:

*   Lopez-Corbalan, R
*   Garcia-Iriepa, C
*   Gramolini, L

Lopez-Corbalan, R

