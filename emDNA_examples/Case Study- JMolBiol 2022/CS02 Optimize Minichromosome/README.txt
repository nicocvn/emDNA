This folder contains files surrounding the production and optimization of a 336-bp minichromosomal circle. Main directory contents including:
- The emDNA_probind command list file: "INPUTFILE"
- The initial condition file: "N336.par"
- The step-parameter file of the 141-bp pathway of the Protein Data Bank file 1KX5: "n141.par"
- The circular coding sequence file associated with N336: "seq_N336_circle.txt"
- A bash script that can be run to execute all emDNA commands for this case study: "cs02.sh"

There are three results directories users can refer to for comparisons.
- "_results-emDNA_probind" contains the three anticipated output files produced once the linear ramping is complete. Note that the iteration count in the "emDNA_pb_stats" file may differ.

- "_results-optimization" contains the optimization log file and optimized step-parameter file for the sequence-dependent optimization.

- "_results-optional-optimization" contains the optimization log file and optimized step-parameter file for the sequence-dependent optimization.

In addition, the "try_nucleosome-sliding" directory contains the code users can run to try the sliding work described at the end of the Case Study section

The contents of this folder are discussed in "emDNA – A Tool for Modeling Protein-decorated DNA Loops and Minicircles at the Base-pair Step Level" found in the 2022 Computational Resources Special Issue of the Journal of Molecular Biology (doi pending).

Please refer to the Case Studies portion of the above article for details regarding this case study.
