#! /bin/bash
# beware: this script is written to be used with the sge gird engine
# if you are running this on your local computer or on a different 
# grid engine the call on line 25 should be changed to `bash` 
# or something else

config_file=$1 # path to the yaml file
output_dir=$2 # where to write the output json files
structure=$3 # the name of the system as provided in the yaml flie

hostname
echo ${config_file}
echo ${output_dir}
echo ${structure}


nr_of_states=$(find ${output_dir}/${structure} -mindepth 1 -maxdepth 1 -type d | wc -l)

for i in $(seq 1 1 ${nr_of_states})
    do 
    for j in $(seq 1 1 ${nr_of_states})
        do 
        echo $i
        echo $j 
        qsub calculate_energies.sh ${config_file} ${output_dir} ${structure} ${i} ${j}
    done
done