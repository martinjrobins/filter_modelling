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
fibres_charges=(0.0001 0.001 0.01 0.1)

echo "blah"
COUNTER=0
for electrostatics_fibre in 1 0
do
    if [ $electrostatics_fibre -eq 1 ] 
    then
        fibres_charges=(0.0001 0.001 0.01 0.1)
    else
        fibres_charges=(0.0)
    fi 
    for fibres_charge in ${fibres_charges[@]}   
    do
        dirname="electro_${electrostatics_fibre}_${fibres_charge}"
        mkdir -p $dirname
        echo "electro = $electrostatics_fibre charge = $fibres_charge"
        echo "doing regular..."
        echo -e "#!/bin/bash'\n'./simulation --fibre_number ${fibre_number} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 0 --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --domain_size ${fibre_number} --base_dir "${dirname}/regular_" &> ${dirname}/run_sim_regular.out" > job_$COUNTER.sh
        COUNTER=$((COUNTER + 1))
        echo "doing hexagon..."
        echo -e "#!/bin/bash'\n'./simulation --fibre_number ${fibre_number} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 1 --domain_size ${fibre_number} --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --base_dir "${dirname}/hexagon_" &> ${dirname}/run_sim_hexagon.out" > job_$COUNTER.sh 
        COUNTER=$((COUNTER + 1))
        for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
        do
            echo "doing ${i}random..."
            echo -e "#!/bin/bash\n./simulation --seed $i --fibre_number ${fibre_number_random} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} --fibre_arrangement 2 --domain_size ${fibre_number} --electrostatics_fibre ${electrostatics_fibre} --fibres_charge ${fibres_charge} --base_dir "${dirname}/${i}random_" &> ${dirname}/run_sim_${i}rand.out" > job_$COUNTER.sh
            COUNTER=$((COUNTER + 1))
        done
    done
done

