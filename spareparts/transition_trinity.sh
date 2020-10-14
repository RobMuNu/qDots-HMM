#!/bin/bash
# Created by Aidan Zecevik
# Last updated by Roberto Munoz 24/08/20
# Make executable with $ chmod 755 runCSSR-allQD-CG2k.sh


# INPUT: 		Path to desired dot file
# OUTPUT: 		MATLAB-readable transition tensor & flat matrix for asymptotics
#
# ASSUMPTIONS:	CG2k-compressed alphabet
#				!!! Almost impossible to tell if alphabet letters have been skipped !!! 

# This grabs the initial probability vector (in a column file) 
# and the transition matrix (for each emitted symbol in a new file) from the CSSR output files 



# Rob's notes
# Note that the nStates can be grabbed from the last transition in the .dot file
# dot will list the FROM states in order, so the last line will contain the
# transition from the larget labelled state. STATE COUNTING STARTS AT ZERO!


# Why reinvent the wheel? Just hack aidans code tbh


transition_grab () {

# Get number of states
nStates=$(tail -2 $1 | head -1 | cut -d " " -f 1)




}

# Run like: transition_grab <argument>


# Read second last line of  file
# tail -2 yourfile | head -1

