#!/bin/bash

# for maps
c0_min=0.01
c0_max=0.5
# for compact
# c0_min=1.0
# c0_max=30.0
nc0=40

nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 10 --filename optimc_nx10.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 15 --filename optimc_nx15.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 20 --filename optimc_nx20.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 25 --filename optimc_nx25.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 30 --filename optimc_nx30.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 35 --filename optimc_nx35.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 40 --filename optimc_nx40.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 45 --filename optimc_nx45.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 50 --filename optimc_nx50.out &
