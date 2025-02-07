#!/usr/bin/env python3

# Author: Lorenzo Gramolini
# GitHub: github.com/lorenzo-gram
# Citation: If you use this code, please cite Gramolini, L


###################################################################################################################
###################################################################################################################
################################### ANALYSIS DOCKING GENERAL  #####################################################
###################################################################################################################
###################################################################################################################
# We select the poses that we want to filter and we remove those who have the same rmsd keeping the one of lower energy inside the same docking.
# We merge the selected poses in a single file, merged.pdbqt
import re
import glob
file_list = glob.glob('*.pdbqt') #glob module is used to generate a list with all the files with .pdbqt extension

pose_inicial=1   #the user defines the first pose to analyse (they're already ordered by energy)
pose_final=10  #the user defines the last pose to analyse (they're already ordered by energy)
pose_inicial=int(pose_inicial)
pose_final=int(pose_final)

for file_name in file_list:    #open the files 1 by 1
    PDBQT=open(file_name,'r')
    pdbqt=[]
    for line in PDBQT:
        splitted=line.split()
        pdbqt.append(splitted)
    
    
##############  ENERGIES  ##########  I want all the energies (affinities) in a list
    energies=[]
    for i in range(len(pdbqt)):
        if len(pdbqt[i])>0 and pdbqt[i][0]=='REMARK' and pdbqt[i][1]=='VINA' and pdbqt[i][2]=='RESULT:':
            energies.append(pdbqt[i][3]) #the energy is in the line starting by "REMARK VINA RESULT:" (1st element after RESULT:)
            
    for i in range(len(energies)):
        energies[i]=float(energies[i])
    
    if pose_inicial==1:
        index_to_delete_for_E = list(range(pose_final, len(energies))) #I create a list with the index of the poses I don't want
    else:                                                              #aka that are not included in [pose_inicial:pose_final] (in Py numeration)
        index_to_delete_for_E = list(range(0, pose_inicial-1)) + list(range(pose_final, len(energies)))

    energies=energies[pose_inicial-1:pose_final] #I keep only the energies of the poses that I want 
    
############  INDECES  ############ These are the indices of the poses I want (in Vina numeration, not Python's that start from 0)
    indices=[]
    for i in range(pose_inicial,pose_final+1):
        indices.append(i)
    
    
#########  X COORDINATE  ######### I extract the X coordinate of the selected poses
    X=[]
    x=[]
    for i in range(len(pdbqt)-1):
        if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
            x.append(pdbqt[i][5])
        if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
            X.append(x)
            x=[]
    
    X=X[pose_inicial-1:pose_final]  #I keep only the coordinate of the selected poses
    
########## Y COORDINATE  #########  Same for y coordinate
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
    
    Y=Y[pose_inicial-1:pose_final]
###########  Z COORDINATE  #########  Same for z coordinate
    Z=[]
    z=[]
    for i in range(len(pdbqt)-1):
        if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
            z.append(pdbqt[i][7])
        if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
            Z.append(z)
            z=[]
    Z=Z[pose_inicial-1:pose_final]
    
##########  ATOM SYMBOL  ######### I extract the atom symbol for all the atoms in the selected poses
    ATOM=[]
    atom=[]
    for i in range(len(pdbqt)-1):
        if len(pdbqt[i])>0 and pdbqt[i][0]=='HETATM':
            atom.append(pdbqt[i][2])
        if len(pdbqt[i+1])>0 and pdbqt[i+1][0]=='ENDMDL':
            ATOM.append(atom)
            atom=[]
    ATOM=ATOM[pose_inicial-1:pose_final]
    
    
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
                sum_of_dist.append(dist_tot_less) #the smallest distance is stored in dist_tot_less and the k loop moves the the 
                                                  #following atom
    
            rmsd=(sum(sum_of_dist)/len(sum_of_dist))**0.5  #this is rmsd'(i,j)
    
            RMSD.append([i+pose_inicial,j+pose_inicial,rmsd]) #RMSD has the indices of the poses in the considered pair and their rmsd'
    
    
    RMSD_final=[]
    for i in range(len(RMSD)):
        for j in range(i,len(RMSD)):
            if RMSD[i][0]==RMSD[j][1] and RMSD[i][1]==RMSD[j][0]: # I'm looking for the pairs (i,j) and (j,i)
                RMSD_final.append([RMSD[i][0],RMSD[i][1],max(RMSD[i][2],RMSD[j][2])])
    #print(RMSD_final)
    
    
################ CREATE NEW PDBQT ONLY WITH THE SELECTED POSES ##########################################
    
    
##########  CRITERION FOR RMSD ##########################
    index_to_delete_for_RMSD=[]
    for i in range(len(RMSD_final)):
        if RMSD_final[i][2]<2:   ####### !!!!!!!!!!     HERE 2 IS USER DEFINED PARAMETER    !!!!!!!! #######################
            if energies[RMSD_final[i][0]-pose_inicial]<energies[RMSD_final[i][1]-pose_inicial]: # take the lower energy pose between 2 poses
                if RMSD_final[i][1]-1 not in index_to_delete_for_RMSD:                          # with same rmsd 
                    index_to_delete_for_RMSD.append(RMSD_final[i][1]-1)
    
            if energies[RMSD_final[i][0]-pose_inicial]>energies[RMSD_final[i][1]-pose_inicial]: 
                if RMSD_final[i][0]-1 not in index_to_delete_for_RMSD:
                    index_to_delete_for_RMSD.append(RMSD_final[i][0]-1)
    
    index_to_delete_for_RMSD.sort()
    
#####  MERGE CRITERIA  #### MERGE THE CRITERIA FOR ENERGY AND RMSD
    delete=list(set(index_to_delete_for_E)|set(index_to_delete_for_RMSD))
    
############ FIND THE INDICES TO DELETE #######################
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
    length=b-a   ### length (in lines) of one pose in the pdbqt
    
    delete_lines=[]
    for i in range(len(delete)):
        for j in range(delete[i]*length+delete[i],(delete[i]+1)*length+delete[i]+1): #the indices of all the lines containing the ith pose
            delete_lines.append(j)
    #print(delete_lines)
    
    
############ WRITE NEW FILE  #############################################
    with open(file_name, 'r') as file:
        # Read all lines from the input file
        lines = file.readlines()
    
    with open('selected_poses_{}'.format(file_name), 'w') as file:
        # Iterate through each line with its index
        for index, line in enumerate(lines):
            # Check if the current index is not in the list of indices to delete
            if index not in delete_lines:
                # Write the new file only with the lines to keep
                file.write(line)
    
######################################################################################################################    
##########################MERGE ALL THE FILES WITH THE SELECTED POSES FOR EACH DOCKING#####################
###############################################################################################################    

# Find all files whose name start with "selected_poses_"
import glob

file_pattern = "selected_poses_*"
files = glob.glob(file_pattern)

# Apri il file di output
with open('merged.pdbqt', 'w') as outfile:
    for filename in files:
        with open(filename, 'r') as infile:
            content = infile.read().rstrip()  # Remove spaces at the end of each line, if any
            outfile.write(content + '\n')  # write each file under the other without blank lines between them


    
    
    
