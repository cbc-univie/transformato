states=$1 # number intst dirs - 1
r=${2:-1}
method=${3:-1}
#blk=10000

mkdir -p ./tmp
rm -f ./tmp/prod_${method}.${r}.output
rm -f ./tmp/raw*dat

for (( i=1; i<=${states}; i++ ))
do
    j=$(($i+1));
    echo "$i -> $j"
    paste intst${i}/prod_${method}.$r.${i}_${i}.dat intst${j}/prod_${method}.$r.${j}_${i}.dat \
	  intst${j}/prod_${method}.$r.${j}_${j}.dat intst${i}/prod_${method}.$r.${i}_${j}.dat > ./tmp/raw.${r}.$i.dat;
    ~/work/festutorial/absolute/newbar/newbar1 -f ./tmp/raw.${r}.$i.dat | tee -a ./tmp/prod_${method}.${r}.output

#   the subblocks make little sense for this systeme 
#    for idx in {1..4}
#    do
#	head -$(($idx*${blk})) tmp/fw.dat | tail -${blk} > tmp/fw.tmp
#	head -$(($idx*${blk})) tmp/bw.dat | tail -${blk} > tmp/bw.tmp
#	luajit ../template_scripts/bar1.lua tmp/fw.tmp tmp/bw.tmp | awk '/dabar/ {print $3}' >> data/part${idx}.${r}.result;
#	echo "------"
#    done
done

echo ""
echo "== Cum. dA =="
echo 
awk '/da_est/ {sum+=$2; i++; printf("%3d%10.3f%10.3f\n",i,$2,sum)}' ./tmp/prod_${method}.${r}.output | tee ./data/prod_${method}.${r}.dat

