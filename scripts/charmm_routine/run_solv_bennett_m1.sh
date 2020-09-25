molecule=$1
idxmax=$2
#state=$2
#rand=$3
method=1 #hardwired
#===========
cd $molecule
for (( idx=1; idx<=${idxmax}; idx++ ))
do
    dir="intst${idx}"
    echo $dir
    cd $dir
    for i in {1..5}
    do
	#echo "charmm mol:$molecule idx:$idx idxmax:$idxmax rand:$i -i get_gasp_energies.inp &> /dev/null"
	charmm mol:$molecule idx:$idx idxmax:$idxmax method:$method rand:$i -i get_energies.inp &> /dev/null &
    done
    wait
    cd ..
done
