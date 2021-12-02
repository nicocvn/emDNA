#!/bin/bash

# ---
# This script will do three tasks (Steps 1, 2, 3) as described below and in the Journal of Molecular Biology article mentioned below.
# To run this script in the same folder or directory as the initial file, type the following to the command line: $ ./cs02.sh
#
# -!-Author note: Based on the build of the user's machine, running this script may take over an hour. Please plan accordingly.
# 
# Script made by Robert Young, 2021
# Script tied to the Journal of Molecular Biology article: emDNA – A Tool for Modeling Protein-decorated DNA Loops and Minicircles at the Base-pair Step Level
# ---

FILENAME=N336_n141

# Step 1: Run the linear ramping function
# --- Author suggestion: Double-check the contents of the INPUTFILE such that all commands required are correct
emDNA_probind INPUTFILE

# --- clean up data files into consistent names
cp emDNA_pb_confs.txt ${FILENAME}_confs.txt
cp emDNA_pb_final.txt ${FILENAME}.txt
cp emDNA_pb_stats.txt ${FILENAME}_stats.txt



# Step 2: Parse the ramping product from a reference frame list to a step-parameter file with sequence information
# --- Prepare the sequence information needed for Step 2
sequencefile='seq_N336_circle.txt'
while read line; do SEQ=${line}; done < ${sequencefile}

emDNA_parser \
--bp-list-input=${FILENAME}.txt \
--bp-list-sequence=${SEQ} \
--get-x3DNA-params>${FILENAME}.par



# Step 3: Optimize the emDNA_produce using a sequence-dependent elastic parameter model
emDNA \
 --x3DNA-bp-step-params-input=${FILENAME}.par \
 --DNA-seqdep-model=Olson1998 \
 --frozen-steps=1:141 \
 --hold-last-bp \
 --output-name=${FILENAME}_seqdep_circle



### OPTIONAL COMMANDS ---

# Option 1: Run free-collection on the emDNA_probind product to see how the linear piece compares with the bound nucleosomal pathway
#emDNA \
# --x3DNA-bp-step-params-input=${FILENAME}.par \
# --DNA-seqdep-model=Olson1998 \
# --frozen-steps=1:141 \
# --free-collection \
# --output-name=${FILENAME}_seqdep_free

# Option 2: Optimized the emDNA_probind product with the homopolymeric IdealDNA elastic model to get a detailed log file of this optimized product
# --- this will be a quick run as emDNA_probind used IdealDNA to optimized
#emDNA \
# --x3DNA-bp-step-params-input=${FILENAME}.par \
# --DNA-seqdep-model=IdealDNA \
# --frozen-steps=1:141 \
# --hold-last-bp \
# --output-name=${FILENAME}_ideal_circle