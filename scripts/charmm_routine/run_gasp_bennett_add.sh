molecule=$1
idxini=$2
idxmax=$3
#state=$2
#rand=$3
#method=$4
#===========
cd $molecule
for (( idx=${idxini}; idx<=${idxmax}; idx++ ))
do
    dir="intst${idx}"
    echo $dir
    cd $dir
    for i in {1..5}
    do
	#echo "charmm mol:$molecule idx:$idx idxmax:$idxmax rand:$i -i get_gasp_energies.inp &> /dev/null"
	charmm mol:$molecule idx:$idx idxmax:$idxmax rand:$i -i get_gasp_energies.inp &> /dev/null &
    done
    wait
    cd ..
done
