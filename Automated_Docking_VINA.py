#!/usr/bin/env python3

# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
# Citation: If you use this code, please cite Lopez-Corbalan, R

import os
import shutil
import subprocess
import time
import glob



# Get the name of the current folder (which will be the parent folder)
parent_folder = os.path.basename(os.getcwd())

# Launch script.py
print("Launching docking_calculus.py")
subprocess.run(["python3", "docking_calculus.py"])
print("docking_calculus.py has finished")

# Search for the only folder in the current directory
folders = [d for d in os.listdir() if os.path.isdir(d)]
if len(folders) != 1:
    raise Exception("Expected to find only one folder, but found {}".format(len(folders)))

folder = folders[0]

# Continue with the rest of the script
print(f"script.py has finished for {folder}")
print(f"Performing the rest of the processes in {folder}")

# Change to the directory of the found folder
os.chdir(folder)

# Copy all necessary scripts from the previous folder
print(f"Copying all scripts to {folder}")
shutil.copy("../add_name.py", ".")
shutil.copy("../analysis_docking_general_auto.py", ".")
shutil.copy("../analysis_docking_particular_auto.py", ".")
shutil.copy("../sort_poses.py", ".")
print(f"All scripts have been successfully copied to {folder}")

# Execute add_name.py
print("add_name.py is about to be launched")
subprocess.run(["python3", "add_name.py"])
print("add_name.py has finished")

# Execute analysis_docking_general_auto.py
print("analysis_docking_general_auto.py is about to be launched")
subprocess.run(["python3", "analysis_docking_general_auto.py"])
print("analysis_docking_general_auto.py has finished")

# Creating final_selection_1
print("final_selection_1 will be created and analysis_docking_particular_auto.py, sort_poses.py will be copied, and the merged.pdbqt file will be moved")
os.mkdir("final_selection_1")
os.chdir("final_selection_1")

# Move necessary files to final_selection_1
shutil.copy("../analysis_docking_particular_auto.py", ".")
shutil.copy("../sort_poses.py", ".")
shutil.move("../merged.pdbqt", ".")
print("Files correctly placed in final_selection_1")

# Launch analysis_docking_particular_auto.py
print("analysis_docking_particular_auto.py is about to be launched")
subprocess.run(["python3", "analysis_docking_particular_auto.py"])
print("analysis_docking_particular_auto.py has finished")

# Create the ultimate_selection folder and move files
print("The ultimate_selection folder will be created, necessary .pdbqt files will be copied, and sort_poses.py will be launched")
os.mkdir("ultimate_selection")
os.chdir("ultimate_selection")

# Move .pdbqt files
for file in glob.glob("../*_*.pdbqt"):
    shutil.move(file, ".")

if os.path.exists("final_selection.pdbqt"):
    os.remove("final_selection.pdbqt")

shutil.move("../sort_poses.py", ".")
print("The ultimate_selection folder has been created and all necessary files have been moved to it")

# Launch sort_poses.py
print("sort_poses.py is about to be launched")
subprocess.run(["python3", "sort_poses.py"])
print("sort_poses.py has finished")

# Completion notification
print(f"Complete process finished in {folder}")
