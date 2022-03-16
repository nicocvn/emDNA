#!/bin/bash

#--- To be run once an initial condition is made.
#--- Make initial configuration parameter of 336-bp with an even-Lk topology using "cs00" Python3 script
python planar_circle_generator.py -l 336 -s seq_zech336.txt -t 10.5 -o N336

#--- Make initial configuration parameter of 195-bp with an even-Lk topology using "cs00" Python3 script
python planar_circle_generator.py -l 195 -s seq_circ195.txt -t 10.5 -o circ195
