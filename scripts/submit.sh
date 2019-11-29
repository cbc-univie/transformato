#! /bin/bash

config_file=$1
working_dir=$2
structure=$3

hostname
echo ${config_file}
echo ${working_dir}
echo ${structure}


nr_of_states=$(find ${working_dir}/${structure} -mindepth 1 -maxdepth 1 -type d | wc -l)

for i in $(seq 1 1 ${nr_of_states})
    do 
    for j in $(seq 1 1 ${nr_of_states})
        do 
        echo $i
        echo $j 
        qsub calculate_energies.sh ${config_file} ${working_dir} ${structure} ${i} ${j}
    done
done