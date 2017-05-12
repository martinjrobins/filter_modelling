#!/bin/bash

# for maps
nbucket_min=10
nbucket_max=20
fr=0.1
nx=10
c0=0.0345

nohup nice ./fmm_eval --nbucket_min ${nbucket_min} --nbucket_max ${nbucket_max} --nx ${nx} --fibre_resolution ${fr} --c0 ${c0} &

