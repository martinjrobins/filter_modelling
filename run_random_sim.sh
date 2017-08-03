#!/bin/bash

# for maps
fr=0.2
c0=0.0835
#fr=0.1
#c0=0.0345
nx=20
fibre_number=10
fibre_number_random=100
electrostatics_fibre=0
fibres_charge=0.1

echo "blah"

for i in 01 02 03 04 05 06 07 08 09 10
do
    nohup nice ./simulation --seed $i --fibre_number ${fibre_number_random} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 2 --domain_size ${fibre_number} --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --base_dir "${i}random_" &> run_sim_${i}rand.out &
done
