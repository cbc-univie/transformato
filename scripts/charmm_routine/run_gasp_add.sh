molecule=$1
#state=$2
#rand=$3
#method=$4
#===========
cd $molecule
for dir in intst1[678]
do
    echo $dir
    cd $dir
    for i in {1..5}
    do
	charmm mol:$molecule rand:$i -i run_gasp_md.inp > run_gasp_md.$i.out&
    done
    wait
    cd ..
done
