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

matrix_N_list=(250 500 750 1250 1500)
process_num_list=(1 4 8 12 16 20 24 28 32)
echo matrix_N: ${matrix_N_list[@]}
echo process_num_list: ${process_num_list[@]}

for matrix_N in ${matrix_N_list[@]}
do
    i=0
    echo "running $(( ${matrix_N} )) matrix_N"
    python3 gen_matrix.py ${matrix_N}
    for process_num in ${process_num_list[@]}
    do
        echo "running $(( ${process_num} )) processes"

        #modify hostfile
        export OMP_NUM_THREADS=${process_num}
        num_list[$i]=${process_num}

        output=$(${program} ${program_args})
        re="elapsed time:\ *([0-9\.]+)"
        if [[ ${output} =~ ${re} ]]
        then
            result[$i]=${BASH_REMATCH[1]}
        else
            echo "Regex match failed, exit"
            exit 0
        fi

        ((i=i+1))
    done
    echo result: ${result[@]}
    echo num_list: ${num_list[@]}
done
