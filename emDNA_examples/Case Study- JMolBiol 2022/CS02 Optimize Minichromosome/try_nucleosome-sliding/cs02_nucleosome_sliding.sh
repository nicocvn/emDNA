#!/bin/bash

# This is an optional method that expands on Case Study 2 whose data is mentioned in the Case Study section of the Journal of Molecular Biology article below.
# In order to run this script, you will need the emDNA_probind product in its original base pair list file format.
# The end result of this script will be to produce 336 unique minicircles where the N336 sequence (see: Figure 2A, JMB article) moves along the DNA pathway, simulating protein sliding along the fragment.
# To run this code: $ ./cs02_nucleosome_sliding.sh
#
# -!-Author note: Based on the build of the user's machine, running this script may take over an hour. Please plan accordingly.
# 
# Script made by Robert Young, 2021
# Script tied to the Journal of Molecular Biology article: emDNA – A Tool for Modeling Protein-decorated DNA Loops and Minicircles at the Base-pair Step Level
# ---

# 1. Make a new directory that will store all of the optimized data
mkdir optimization_data

# 2. Prepare the variables for the code and upload the sequence string
FILENAME=N336_n141
sequencefile='seq_zech336.txt'
while read line; do SEQUENCE=${line}; done < ${sequencefile}

# Optimization code block
for ((i=1; i<=${#SEQUENCE}; i++));
do

 # Since 336 minicircles will be made, this line is used to ensure the same filename lengths
 printf -v SLIDEID "%03d" ${i}
 
 # Make sure you are running the correct ID number
 echo ${SLIDEID}
 
 # This is needed to "slide" the sequence by moving 'i' number of characters from the 3` end to the 5` end
 SEQ=${SEQUENCE:i}${SEQUENCE::i}${SEQUENCE::1}

 # This will generate a new step-parameter file with the newly-slid sequence
 emDNA_parser \
 --bp-list-input=${FILENAME}.txt \
 --bp-list-sequence=${SEQ} \
 --get-x3DNA-params>${FILENAME}_seq${SLIDEID}.par

 # Optimization of the minichromosome with the new sequence
 emDNA \
 --x3DNA-bp-step-params-input=${FILENAME}_seq${SLIDEID}.par \
 --DNA-seqdep-model=Olson1998 \
 --frozen-steps=1:141 \
 --hold-last-bp \
 --output-name=${FILENAME}_seq${SLIDEID}_seqdep_circle
 
 # Move the data to the optimization_data directory
 mv *_circle_opt* optimization_data
 
 # -OPTIONAL- Delete the initial file used for optimization_data
 # rm ${FILENAME}_seq${SLIDEID}.par

done


# TROUBLESHOOTING: --free-collection OPTION
#mkdir free_collection_data
#FILENAME=N336_n141
#sequencefile='seq_zech336.txt'
#while read line; do SEQUENCE=${line}; done < ${sequencefile}
#for ((i=1; i<=${#SEQUENCE}; i++));
#do
# printf -v SLIDEID "%03d" ${i}
# echo ${SLIDEID}
# SEQ=${SEQUENCE:i}${SEQUENCE::i}${SEQUENCE::1}
# emDNA_parser \
# --bp-list-input=${FILENAME}.txt \
# --bp-list-sequence=${SEQ} \
# --get-x3DNA-params>${FILENAME}_seq${SLIDEID}.par
# emDNA \
# --x3DNA-bp-step-params-input=${FILENAME}_seq${SLIDEID}.par \
# --DNA-seqdep-model=Olson1998 \
# --frozen-steps=1:141 \
# --free-collection \
# --output-name=${FILENAME}_seq${SLIDEID}_seqdep_free
# mv *_free_* free_collection_data
#done