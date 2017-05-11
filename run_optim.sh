#!/bin/bash

# for maps
c0_min=0.01
c0_max=0.5
# for compact
# c0_min=1.0
# c0_max=30.0
nc0=40

nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 10 --fibre_resolution 0.4 --filename optimc_nx10_f0_4.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 10 --fibre_resolution 0.2 --filename optimc_nx10_f0_2.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 10 --fibre_resolution 0.1 --filename optimc_nx10_f0_1.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 10 --fibre_resolution 0.05 --filename optimc_nx10_f0_05.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 20 --fibre_resolution 0.8 --filename optimc_nx20_f0_8.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 20 --fibre_resolution 0.4 --filename optimc_nx20_f0_4.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 20 --fibre_resolution 0.2 --filename optimc_nx20_f0_2.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 20 --fibre_resolution 0.1 --filename optimc_nx20_f0_1.out &

