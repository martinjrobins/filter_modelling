#!/bin/bash

# for maps
fr=0.2
c0=0.0835
#fr=0.1
#c0=0.0345
nx=20
fibre_number=10
fibre_number_random=100
electrostatics_fibre=1
fibres_charge=0.1

echo "blah"
nohup nice ./simulation --fibre_number ${fibre_number} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 0 --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --domain_size ${fibre_number} --base_dir "regular_" &> run_sim_regular.out &
nohup nice ./simulation --fibre_number ${fibre_number} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 1 --domain_size ${fibre_number} --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --base_dir "hexagon_" &> run_sim_hexagon.out &
nohup nice ./simulation --fibre_number ${fibre_number_random} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 2 --domain_size ${fibre_number} --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --base_dir "random_" &> run_sim_random.out &

