molecule=$1
state=$2
#rand=$3
#method=$4
#===========
cd $molecule/intst${state}
pwd
for method in {1..3}
do
    for i in {1..5}
    do
	echo "Working on method $method and random seed $i .."
	charmm mol:$molecule rand:$i method:${method} -i run_cpt_watmd.inp > run_cpt_watmd_m${method}.$i.out
	charmm mol:$molecule rand:$i method:${method} -i run_cpt_watmd_prod1.inp > run_cpt_watmd_prod1_m${method}.$i.out
    done
done

