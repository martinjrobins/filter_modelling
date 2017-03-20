#!/bin/bash

c0_min=0.01
c0_max=0.5
nc0=40

nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 20 --filename optimc_nx20.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 25 --filename optimc_nx25.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 30 --filename optimc_nx30.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 35 --filename optimc_nx35.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 40 --filename optimc_nx40.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 45 --filename optimc_nx45.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 50 --filename optimc_nx50.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 55 --filename optimc_nx55.out &
nohup nice ./optimc --c0_min ${c0_min} --c0_max ${c0_max} --nc0 ${nc0} --nx 60 --filename optimc_nx60.out &
