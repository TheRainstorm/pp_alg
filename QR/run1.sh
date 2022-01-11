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

process_num_list=(1 4 8 12 16 20 24 28 32)
echo process_num_list: ${process_num_list[@]}

i=0
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
        echo $result[$i]
    else
        echo "Regex match failed, exit"
        exit 0
    fi

    ((i=i+1))
done
echo result: ${result[@]}
echo num_list: ${num_list[@]}

