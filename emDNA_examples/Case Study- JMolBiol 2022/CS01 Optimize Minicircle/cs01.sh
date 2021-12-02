#!/bin/bash

# ---
# This script will run emDNA on the N336.par initial configuration using the homopolymeric and the sequence dependent elastic models.
# To run this script in the same folder or directory as the initial file, type the following to the command line:
# $ ./cs01.sh
#
# Script made by Robert Young, 2021
# Script tied to the Journal of Molecular Biology article: emDNA – A Tool for Modeling Protein-decorated DNA Loops and Minicircles at the Base-pair Step Level
# ---



FILENAME=N336

emDNA \
 --x3DNA-bp-step-params-input=${FILENAME}.par \
 --DNA-seqdep-model=IdealDNA \
 --hold-last-bp \
 --output-name=${FILENAME}_ideal_circle

emDNA \
 --x3DNA-bp-step-params-input=${FILENAME}.par \
 --DNA-seqdep-model=Olson1998 \
 --hold-last-bp \
 --output-name=${FILENAME}_seqdep_circle



### METHOD BELOW IS FOR TROUBLESHOOTING: RUN --free-collection OPTION
#emDNA \
# --x3DNA-bp-step-params-input=${FILENAME}.par \
# --DNA-seqdep-model=IdealDNA \
# --free-collection \
# --output-name=${FILENAME}_ideal_free

#emDNA \
# --x3DNA-bp-step-params-input=${FILENAME}.par \
# --DNA-seqdep-model=Olson1998 \
# --free-collection \
# --output-name=${FILENAME}_seqdep_free


