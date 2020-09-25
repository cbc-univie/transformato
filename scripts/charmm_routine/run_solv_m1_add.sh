molecule=$1
#state=$2
method=1
#rand=$3
#method=$4
#===========
cd $molecule

for state in intst1[678]
do
    cd $state
    for i in {1..5}
    do
	echo "Working on method $method and random seed $i, state $state .."
	OMP_NUM_THREADS=8 charmm mol:$molecule rand:$i method:${method} -i run_cpt_watmd.inp > run_cpt_watmd_m${method}.$i.out
	OMP_NUM_THREADS=8 charmm mol:$molecule rand:$i method:${method} -i run_cpt_watmd_prod1.inp > run_cpt_watmd_prod1_m${method}.$i.out
    done
    cd ..
done

