#!/bin/bash

# for maps
fr=0.1
c0=0.0345
nx=20
fibre_number=10

nohup nice ./simulation --fibre_number ${fibre_number} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 0 --base_dir "regular_" &> run_sim_regular.out &
nohup nice ./simulation --fibre_number ${fibre_number} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 1 --base_dir "hexagon_" &> run_sim_hexagon.out &

