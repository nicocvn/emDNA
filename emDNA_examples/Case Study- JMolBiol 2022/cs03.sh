#!/bin/bash

#-!- TO RUN THIS CODE, USERS MUST HAVE THE FOLLOWING EITHER IN THEIR PATH OR IN CURRENT DIRECTORY
#--- circ195.par (made using "cs00" Python3 script)
#--- n141.par (included in supplied data)
#--- seq_circ195.txt (includes bp148-bp6 in sequence in JMB 2022 paper in Figure 2A OR included file in supplied data)


# ----------------------------------------------------------------------------------------

# --- Generate loop ends from the n141.par initial condition
# -!- Optional organization code
mkdir initial_loop_ends

# --- extract the first and last base-pair frames
emDNA_parser \
 --x3DNA-bp-step-params-input=n141.par \
 --get-bp-list>n141.txt
 
echo $(head -n 1 n141.txt) >> step.txt
echo $(tail -n 1 n141.txt) >> step.txt

emDNA_parser \
 --bp-list-input=step.txt \
 --bp-list-sequence='AA' \
 --get-x3DNA-params>step_std.par
 
rm n141.txt step.txt

# -!- Take the step parameter and replace the RISE parameter with values ranging
# -!- ... from 0-Ang to -500-Ang in steps of -50-Ang

for ((i=0; i<=10; i++));
do
 RISE=$(printf %.4f "$((i*-50))")
 lastpar=$(tail -n 1 step_std.par)
 read -a strarr <<< $lastpar
 sed "s/${strarr[9]}/$RISE/" step_std.par > step_${i}.par
 unset RISE
 unset lastpar
 mv step_${i}.par initial_loop_ends
done
mv step_std.par initial_loop_ends

# ----------------------------------------------------------------------------------------
# --- Generate loop ends from the n141.par initial condition
# -!- Optional organization code
mkdir initial_data
mkdir probind_data
mkdir optimized_data

# Use linear ramping on the various step parameters produced at various RISE values and
# ... re-optimizing with IdealDNA to collect detailed energetic values
for ((i=0; i<=10; i++));
do
 
 # -*- Note to users: the i value is connected to the rise value, 
 #     e.g. step_3 -> rise = -150.0000
 #     e.g. step_5 -> rise = -250.0000

 MODEL=step_${i}
 
 cp initial_loop_ends/${MODEL}.par .
 
 > input_${MODEL}.txt
 echo 'collection_type=EEDR
 base_collection=x3DNAparams:circ195.par
 protein_collection=x3DNAparams:'${MODEL}'.par
 protein_binding_sites=1
 base_bound_domains=1:2
 binding_ramp_sampling=10
 base_seqdep_model=IdealDNA' >input_${MODEL}.txt

 emDNA_probind input_${MODEL}.txt

 mv emDNA_pb_confs.txt circ195_${MODEL}_confs.txt
 mv emDNA_pb_final.txt circ195_${MODEL}_final.txt
 mv emDNA_pb_stats.txt circ195_${MODEL}_stats.txt

 sed -i '${/^$/d}' circ195_${MODEL}.txt
 
 while read line; do SEQ=${line}; done < 'seq_circ195.txt' 
 CIRCSEQ=${SEQ}${SEQ:0:1}
 
 emDNA_parser \
 --bp-list-input=circ195_${MODEL}_final.txt \
 --bp-list-sequence=${CIRCSEQ} \
 --get-x3DNA-params>circ195_${MODEL}.par

 emDNA \
 --x3DNA-bp-step-params-input=circ195_${MODEL}.par \
 --DNA-seqdep-model=IdealDNA \
 --frozen-steps=1:2 \
 --hold-last-bp \
 --output-name=circ195_${MODEL}

## -!- Optional organization code
 mv ${MODEL}_confs.txt ${MODEL}_final.txt ${MODEL}_stats.txt probind_data
 mv ${MODEL}.log ${MODEL}_opt.txt optimized_data
 mv ${MODEL}.par input_${MODEL}.txt initial_data
 
done
# -!- Optional organization code
mv circ195.par n141.par initial_data
# ----------------------------------------------------------------------------------------
