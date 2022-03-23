#!/bin/bash

#--- To be run once an initial condition is made.
#--- Make initial configuration parameter of 336-bp using "cs00" Python3 script

# -!- Optional organization code
mkdir optimized_data
mkdir initial_data

# --- IdealDNA ---------------------------------------------------------------------------
emDNA \
 --x3DNA-bp-step-params-input=N336.par \
 --DNA-seqdep-model=IdealDNA \
 --hold-last-bp \
 --output-name=N336_ideal_circle

emDNA_topology \
 --x3DNA-bp-step-params-input=N336_ideal_circle_opt.txt>topo_N336_ideal_circle_opt.txt \
 --virtual-last-bp

# --- Olson1998 sequence-dependent forcefield --------------------------------------------
emDNA \
 --x3DNA-bp-step-params-input=N336.par \
 --DNA-seqdep-model=Olson1998 \
 --hold-last-bp \
 --output-name=N336_seqdep_circle

emDNA_topology \
 --x3DNA-bp-step-params-input=N336_seqdep_circle_opt.txt>topo_N336_seqdep_circle_opt.txt \
 --virtual-last-bp

# ----------------------------------------------------------------------------------------
# -!- Optional organization code
mv N336_ideal_circle.log N336_ideal_circle_opt.txt topo_N336_ideal_circle_opt.txt optimized_data
mv N336_seqdep_circle.log N336_seqdep_circle_opt.txt topo_N336_seqdep_circle_opt.txt optimized_data
mv N336.par initial_data
# ----------------------------------------------------------------------------------------

