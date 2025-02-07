#!/usr/bin/env python3

# Author: Lorenzo Gramolini
# GitHub: github.com/lorenzo-gram
# Citation: If you use this code, please cite Gramolini, L


###################################################################################################################
###################################################################################################################
################################### ANALYSIS DOCKING PARTICULAR  #####################################################
###################################################################################################################
###################################################################################################################
# We delete the poses with the same rmsd keeping the one of lower energy in merged.pdbqt
# We create a new file called final_selection.pdbqt with all the selected poses
# We create a .pdbqt file for each pose in final_selection.pdbqt (deleting the lines that cause problems with visualization)
import re
import glob
file_list = glob.glob('*.pdbqt') #glob module is used to generate a list with all the files with .pdbqt extension

PDBQT=open(file_list[0],'r')  #open the files 1 by 1
pdbqt=[]
for line in PDBQT:
    splitted=line.split()
    pdbqt.append(splitted)


##############  ENERGIES  ##########  I want all the energies (affinities) in a list
energies=[]
for i in range(len(pdbqt)):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='REMARK' and pdbqt[i][1]=='VINA' and pdbqt[i][2]=='RESULT:':
        energies.append(pdbqt[i][3])  #the energy is in the line starting by "REMARK VINA RESULT:" (1st element after RESULT:)
        
for i in range(len(energies)):
    energies[i]=float(energies[i])  



######### X COORDINATE ##########################
X=[]
x=[]
for i in range(len(pdbqt)-1):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
        x.append(pdbqt[i][5])
    if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
        X.append(x)
        x=[]


########## Y COORDINATE ############################
Y=[]
y=[]
for i in range(len(pdbqt)-1):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
        if '-' in pdbqt[i][6][1:]:
            match = re.match(r'(-?\d+\.\d+)(-?\d+\.\d+)', pdbqt[i][6])
            if match:
                y.append(match.group(1))
                pdbqt[i][7] = match.group(2)
        else:
            y.append(pdbqt[i][6])
    if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
        Y.append(y)
        y=[]


######### Z COORDINATE  ###########################
Z=[]
z=[]
for i in range(len(pdbqt)-1):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
        z.append(pdbqt[i][7])
    if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
        Z.append(z)
        z=[]

#########  ATOM SYMBOL  ###########################
ATOM=[]
atom=[]
for i in range(len(pdbqt)-1):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
        atom.append(pdbqt[i][2])
    if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
        ATOM.append(atom)
        atom=[]


##########  RMSD lb  ########### I calculate the RMSD lb
RMSD=[]
for i in range(len(ATOM)):
    for j in range(len(ATOM)):
#all the pairwise interaction with permutational order (since rmsd lb = max(rmsd'(1,2),rmsd'(2,1)))
        sum_of_dist=[]  #this list stores the squared distances for all the atom pairs for the current pair of poses
        for k in range(len(ATOM[i])):
            dist_tot_less=10000000  # random number that is always higher than the others
            for t in range(len(ATOM[j])):
                if ATOM[i][k]==ATOM[j][t]:  #if they have the same atom symbol
                #all the combination of atoms
                    dist_x=float(X[i][k])-float(X[j][t])
                    dist_y=float(Y[i][k])-float(Y[j][t])
                    dist_z=float(Z[i][k])-float(Z[j][t])
                    dist_tot=dist_x**2+dist_y**2+dist_z**2
                    if dist_tot<dist_tot_less:
                        dist_tot_less=dist_tot  # it tries the distances between all the atoms in one poses and all the atoms
                                                # in the other pose (with same atom type) and select the smaller distance
            sum_of_dist.append(dist_tot_less)   #the smallest distance is stored in dist_tot_less and the k loop moves the the
                                                #following atom

        rmsd=(sum(sum_of_dist)/len(sum_of_dist))**0.5  #this is rmsd'(i,j)

        RMSD.append([i+1,j+1,rmsd])  #RMSD has the indices of the poses in the considered pair and their rmsd'


RMSD_final=[]
for i in range(len(RMSD)):
    for j in range(i,len(RMSD)):
        if RMSD[i][0]==RMSD[j][1] and RMSD[i][1]==RMSD[j][0]:  # I'm looking for the pairs (i,j) and (j,i)
            RMSD_final.append([RMSD[i][0],RMSD[i][1],max(RMSD[i][2],RMSD[j][2])])
#print(RMSD_final)


################ CREATE NEW PDBQT WITH THE FINAL SELECTION ###################################################


########## CRITERION FOR RMSD ######################################
index_to_delete_for_RMSD=[]
for i in range(len(RMSD_final)):
    if RMSD_final[i][2]<2:   ####### !!!!!!!!!!     HERE 2 IS USER DEFINED PARAMETER    !!!!!!!! #######################
        if energies[RMSD_final[i][0]-1]<energies[RMSD_final[i][1]-1]: # take the lower energy pose between 2 poses
            if RMSD_final[i][1]-1 not in index_to_delete_for_RMSD:    # with same rmsd
                index_to_delete_for_RMSD.append(RMSD_final[i][1]-1)

        if energies[RMSD_final[i][0]-1]>energies[RMSD_final[i][1]-1]:
            if RMSD_final[i][0]-1 not in index_to_delete_for_RMSD:
                index_to_delete_for_RMSD.append(RMSD_final[i][0]-1)

index_to_delete_for_RMSD.sort()
delete=index_to_delete_for_RMSD

############FIND THE INDICES TO DELETE ########################################
a=0
b=0
for i in range(len(pdbqt)):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='MODEL':
        a=i
        break
for i in range(len(pdbqt)):
    if len(pdbqt[i])>0 and pdbqt[i][0]=='ENDMDL':
        b=i
        break
length=b-a

delete_lines=[]
for i in range(len(delete)):
    for j in range(delete[i]*length+delete[i],(delete[i]+1)*length+delete[i]+1):
        delete_lines.append(j)
#print(delete_lines)


########### WRITE NEW FILE #######################################
with open(file_list[0], 'r') as file: # There has to be only one .pdbqt file in the folder where I execute (the merged.pdbqt)
        # Read all lines from the input file
        lines = file.readlines()
    
with open('final_selection.pdbqt', 'w') as file:
    # Iterate through each line with its index
    for index, line in enumerate(lines):
        # Check if the current index is not in the list of indices to delete
        if index not in delete_lines:
            # Write the line to the output file
            file.write(line)

############## REMOVE LINES CONTAINING PROBLEMATIC WORDS FOR VISUALIZATION (BRANCH ETC...) #####################
with open('final_selection.pdbqt', 'r') as file:
    lines = file.readlines()

with open('final_selection.pdbqt', 'w') as file:
    for line in lines:
        if not any(word in line for word in ["ENDROOT", "TORSDOF", "BRANCH", "ENDBRANCH"]):
            file.write(line)

####### SPLIT final_selection.pdbqt INTO A FILE FOR EACH SELECTED POSE #############################

# Open the original file for reading
with open('final_selection.pdbqt', 'r') as file:
    lines = file.readlines()

block = []  # To store the lines of the current block
filename = None  # Placeholder for the filename
file_counter = 0  # To avoid name collisions if 'MODEL' line is missing words

# Iterate through each line of the original file
for line in lines:
    # If the line starts with 'MODEL', begin a new block
    if line.startswith('MODEL'):
        # If there's already a block being built, write it to a file
        if block:
            # Write the previous block to a file
            with open(f'{filename}', 'w') as subfile:
                subfile.writelines(block)
            block = []  # Clear the block for the new one
        
        # Extract the filename from the line (assuming words after 'MODEL' form the name)
        filename_parts = line.split()[1:]  # Extract everything after 'MODEL'
        if filename_parts:
            filename = "_".join(filename_parts)  # Join words with underscores for filename, ITS NAME WILL CONTAIN INFO ON WHICH DOCJING
        else:                                                                              # THE POSE IS FROM
            # Fallback if no extra words after 'MODEL' (to avoid empty filenames)
            file_counter += 1
            filename = f'model_{file_counter}'
        
        # Start a new block with the 'MODEL' line
        block.append(line)
    
    # If the line starts with 'ENDMDL', close the current block
    elif line.startswith('ENDMDL'):
        block.append(line)
        # Write the block to a file since it has finished
        with open(f'{filename}', 'w') as subfile:
            subfile.writelines(block)
        block = []  # Clear the block for the next one
    
    # Otherwise, just keep adding lines to the current block
    else:
        if block:
            block.append(line)

# Edge case: if the file ends without 'ENDMDL', write the remaining block (if any)
if block:
    with open(f'{filename}', 'w') as subfile:
        subfile.writelines(block)




