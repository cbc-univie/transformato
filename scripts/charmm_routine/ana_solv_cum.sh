states=$1 # number intst dirs - 1
method=${2:-1}
#blk=10000

mkdir -p ./tmp
rm -f ./tmp/prod_${method}.cum.output
rm -f ./tmp/raw*dat

for (( i=1; i<=${states}; i++ ))
do
    j=$(($i+1));
    echo "$i -> $j"
    cat intst${i}/prod_${method}.?.${i}_${i}.dat > intst${i}/prod_${method}.cum.${i}_${i}.dat
    cat intst${i}/prod_${method}.?.${i}_${j}.dat > intst${i}/prod_${method}.cum.${i}_${j}.dat
    cat intst${j}/prod_${method}.?.${j}_${i}.dat > intst${j}/prod_${method}.cum.${j}_${i}.dat
    cat intst${j}/prod_${method}.?.${j}_${j}.dat > intst${j}/prod_${method}.cum.${j}_${j}.dat
    paste intst${i}/prod_${method}.cum.${i}_${i}.dat intst${j}/prod_${method}.cum.${j}_${i}.dat \
	  intst${j}/prod_${method}.cum.${j}_${j}.dat intst${i}/prod_${method}.cum.${i}_${j}.dat > ./tmp/raw.cum.$i.dat;
    ~/work/festutorial/absolute/newbar/newbar1 -f ./tmp/raw.cum.$i.dat | tee -a ./tmp/prod_${method}.cum.output

#   the subblocks make little sense for this systeme 
#    for idx in {1..4}
#    do
#	head -$(($idx*${blk})) tmp/fw.dat | tail -${blk} > tmp/fw.tmp
#	head -$(($idx*${blk})) tmp/bw.dat | tail -${blk} > tmp/bw.tmp
#	luajit ../template_scripts/bar1.lua tmp/fw.tmp tmp/bw.tmp | awk '/dabar/ {print $3}' >> data/part${idx}.cum.result;
#	echo "------"
#    done
done

echo ""
echo "== Cum. dA =="
echo 
awk '/da_est/ {sum+=$2; i++; printf("%3d%10.3f%10.3f\n",i,$2,sum)}' ./tmp/prod_${method}.cum.output | tee ./data/prod_${method}.cum.dat

