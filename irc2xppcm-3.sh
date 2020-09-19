#!/usr/bin/env bash

#--------------------------------------------------------------------
# Bo Chen 2020-08-16
# Generates Gaussian input files for XP-PCM calculations from existing
# IRC output file
# Usage: irc2xppcm.sh irc.[log|out]
#--------------------------------------------------------------------

# Get the file name without the extension from the first augument
input=$1
filename="${input%.*}"

# Get the number of atoms from the first coordinate block in the IRC file
numberOfAtoms=$(sed -n '1,/Redundant internal/d;/Recover connectivity/q;p' $input | awk 'END {print NR}')

# Pass the number of atoms and filename as variables to awk.
# The awk code processes the whole IRC file and prints multiple
# .com Gaussin input files for XP-PCM calculation
awk -v nOA="$numberOfAtoms" -v fname="$filename" '

# Some variable initiation 
# start (0 or 1) flags the beginning of a coordinate block
# struct counts the structures in the IRC file
# rxnCoord stores the reaction coordinate of each structure; rc is the 
# index. The first structure is the TS with rxnCoord[1]=0.
# b is a random, large number and its value will be reassigned later
BEGIN { start=0; struct=1; rc=1; rxnCoord[rc]=0; b=10000 }

# Get and store nproc, memory, method & basis set, charge & multiplicity
# Current implementation requires that the keywords line in previous 
# IRC calculation must start with "#p method basis-set ...".
# The first three fields in the keywords line are stored.
/%nproc/ {nproc=$1; next}
/%mem/ {mem=$1; next}
/#p/ {method=$1 " " $2 " " $3; next}
/Multiplicity/ {chargeAndMulti=$3 " " $6; next}

# In the IRC output, the coordinates of TS are printed out first, following 
# the "Redundant internal" line, with dilimiter ",".
# Setting start=1 means the next line is the begining of a coordinate block 
# next ignores the rest of patterns and move to the next record/line
/Redundant internal/ { FS=","; start=1; next; }

# In the IRC output, the coordinates of other structures are printed in the 
# "CURRENT STRUCTURE" block, 6 lines after this line.
# Assign current line number to b, and set start=1 at line b+5.
/CURRENT STRUCTURE/ { b=NR; next; }
NR==b+5 { start=1; next; }

# Assign the coordinate lines as strings to a 2D array geom[struct][atom], 
# the first index being the structure number and the second index being 
# the atom number in the structure
start==1 { 
    # The coordinates of TS (struct=1) are printed 
    # differently in the IRC file, so diferent fields 
    # are extracted and stored
    if(struct==1) {
        geom[struct][++atom]=$1 " " $3 " " $4 " " $5;
    } else {
        geom[struct][++atom]=$2 " " $3 " " $4 " " $5; 
    }
    # When the index atom reaches the number of atoms, meaning 
    # all coordinate lines of one structure have been extracted,
    # reset start=0 so that no new line will be extracted.
    # Increment index struct for the next structure and reset 
    # the atom index to 0.
    # Reset the field seperator/delimiter FS to space, which is
    # what is used for coordinates of other structures other than TS.
    if(atom==nOA) { start=0; struct++; atom=0; FS=" "; next; } 
}

# The following three lines deal with reaction coordinates.
# First, use a varable (sign) with value 1 or -1 to indicate 
# either the forward or reverse direction of the IRC path.
/Point Number  1 in FORWARD/ { sign=1; next; }
/Point Number  1 in REVERSE/ { sign=-1; next; }
# Get the reaction coordinates and store them in an array rxnCoord[]. 
# The first structure is the TS with rxnCoord[1]=0 (assigned in the 
# BEGIN section). Structures on the reverse direction have negative 
# reaction coordinate
/NET REACTION COORDINATE/ { rxnCoord[++rc]=$9*sign; next }

# Match either pattern to get the energies of the structures
# The 2nd pattern for TS and the 1st for all other structures
/SCF Done/ || /Energy From Chk/ { energy[++en]=$5; next }

# Loop through the structure numbers to print out .com file for each structure
END { 
    for (i=1;i<struct;i++) {
        # Constrcut an array with rxnCoord as indexes and stores simple 
	# numeric value from 1 to struct}
	array[rxnCoord[i]]=i 
	# Use asorti function to sort the index of array[], generating a
	# new sortedArray[] with simple numeric index that stores the 
	# sorted rxnCoord values. The third argument tells the function how 
	# to sort, in this case index_numeric_ascend
	asorti(array,sortedArray,"@ind_num_asc")
    }
    # Print output files
    for (i=1;i<struct;i++) {
	# Print cpu, memory, kewords lines to a .com file
	printf("%s\n%s\n%s\n\n", nproc, mem, method) > fname"-"i".com"
	# Print the title line containing rxnCoord and energy of the structure.
	# sortedArray[i] gives the sorted rxnCoord value;
	# array[sortedArray[i]] gives the numeric index according to the 
	# sorted order;
	# energy[ array[sortedArray[i]] ] gives the energy according to 
	# the sorted order
	printf("rxnCoord: %.2f energy: %.8f\n\n", 
	    rxnCoord[array[sortedArray[i]]], 
	    energy[array[sortedArray[i]]]) >> fname"-"i".com"
	# Print charge and multiplicity
	printf("%s\n", chargeAndMulti) >> fname"-"i".com"
	# A second loop to print coordinates, adding a leading space for each
	# coordinate line.
	for(j=1;j<=nOA;j++) { 
	    printf(" "geom[array[sortedArray[i]]][j])"\n" >> fname"-"i".com" 
        }
        # Print an additional black line at the end of file
        printf("\n") >> fname"-"i".com"
    } 
}' $input 



