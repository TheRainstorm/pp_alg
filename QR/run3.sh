#!/bin/bash
if [ $# -lt 1 ]
then
    echo "Usage: $0 <program> [arg]..."
    exit 0
else
    program=$1
    program_args="${@:2}"
    echo "program: ${program}"
    echo "args:${program_args}"
fi

matrix_N_list=(250 500 750 1000 1250 1500)
echo matrix_N: ${matrix_N_list[@]}

i=0
for matrix_N in ${matrix_N_list[@]}
do
    echo "running $(( ${matrix_N} )) matrix_N"
    num_list[$i]=${matrix_N}

    sed -i 's/#define N.*$/#define N    '"${matrix_N}/g" serial_QR.c
    ../0c gcc serial_QR.c

    python3 gen_matrix.py ${matrix_N}

    output=$(${program} ${program_args})
    
    re="elapsed time\ *[:=]\ *([0-9\.]+)"
    if [[ ${output} =~ ${re} ]]
    then
        result[$i]=${BASH_REMATCH[1]}
    else
        echo "Regex match failed, exit"
        exit 0
    fi

    ((i=i+1))
    echo result: ${result[@]}
    echo num_list: ${num_list[@]}
done
