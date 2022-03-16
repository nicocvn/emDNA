#!/bin/bash

#-!- TO RUN THIS CODE, USERS MUST HAVE THE FOLLOWING EITHER IN THEIR PATH OR IN CURRENT DIRECTORY
#--- N336.par (made using "cs00" Python3 script)
#--- n141.par (included in supplied data)
#--- seq_N336.txt (includes sequence in JMB 2022 paper in Figure 2A OR included file in supplied data)


# -!- Optional organization code
mkdir initial_data
mkdir probind_data
mkdir optimized_data

# ----------------------------------------------------------------------------------------

> INPUTFILE.txt
echo 'collection_type=EEDR
base_collection=x3DNAparams:N336.par
protein_collection=x3DNAparams:n141.par
protein_binding_sites=1
base_bound_domains=1:141
binding_ramp_sampling=100
base_seqdep_model=IdealDNA' >INPUTFILE.txt

emDNA_probind INPUTFILE.txt

sed -i '${/^$/d}' emDNA_pb_final.txt

while read line; do SEQ=${line}; done < 'seq_N336.txt'

emDNA_parser \
 --bp-list-input=emDNA_pb_final.txt \
 --bp-list-sequence=${SEQ}${SEQ:0:1} \
 --get-x3DNA-params>N336_n141.par
unset SEQ

# --- IdealDNA ---------------------------------------------------------------------------

emDNA \
 --x3DNA-bp-step-params-input=N336_n141.par \
 --DNA-seqdep-model=IdealDNA \
 --frozen-steps=1:141 \
 --hold-last-bp \
 --output-name=N336_n141_ideal_circle
 
emDNA_topology \
 --x3DNA-bp-step-params-input=N336_n141_ideal_circle_opt.txt>topo_N336_n141_ideal_circle_opt.txt \
 --virtual-last-bp

# --- Olson1998 sequence-dependent forcefield --------------------------------------------
emDNA \
--x3DNA-bp-step-params-input=N336_n141.par \
--DNA-seqdep-model=Olson1998 \
--frozen-steps=1:141 \
--hold-last-bp \
--output-name=N336_n141_seqdep_circle

emDNA_topology \
--x3DNA-bp-step-params-input=N336_n141_seqdep_circle_opt.txt>topo_N336_n141_seqdep_circle_opt.txt \
--virtual-last-bp

# ----------------------------------------------------------------------------------------
# -!- Optional organization code
mv emDNA_pb_confs.txt emDNA_pb_final.txt emDNA_pb_stats.txt probind_data
mv N336_n141_ideal_circle.log N336_n141_ideal_circle_opt.txt topo_N336_n141_ideal_circle_opt.txt optimized_data
mv N336_n141_seqdep_circle.log N336_n141_seqdep_circle_opt.txt topo_N336_n141_seqdep_circle_opt.txt optimized_data
mv N336_n141.par INPUTFILE.txt initial_data

mv N336.par n141.par seq_zech336.txt initial_data
# ----------------------------------------------------------------------------------------

