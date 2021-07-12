#!/bin/bash

SIZE=$1
NSTEPS=16

TTSIZE=`expr ${SIZE} \* ${SIZE} \* ${SIZE}`

./generate ${TTSIZE} stencil_data_${SIZE}x${SIZE}x${SIZE} &&
./stencil_7 stencil_data_${SIZE}x${SIZE}x${SIZE} ${SIZE} ${SIZE} ${SIZE} ${NSTEPS} stencil_answer_7_${SIZE}x${SIZE}x${SIZE}_${NSTEPS}steps &&
./stencil_27 stencil_data_${SIZE}x${SIZE}x${SIZE} ${SIZE} ${SIZE} ${SIZE} ${NSTEPS} stencil_answer_27_${SIZE}x${SIZE}x${SIZE}_${NSTEPS}steps
