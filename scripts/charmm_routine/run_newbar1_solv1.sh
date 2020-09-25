r=${1:-1}
blk=10000

mkdir -p ./tmp
rm -f ./tmp/full.${r}.output
#touch tmp/full.result tmp/part1.result tmp/part2.result tmp/part3.result tmp/part4.result #part5.result


for i in {1..3};
do
    j=$(($i+1));
    echo "$i -> $j"
    paste intst${i}/alldat.1.${i}_${i}.$r.dat intst${j}/alldat.1.${j}_${i}.$r.dat \
	  intst${j}/alldat.1.${j}_${j}.$r.dat intst${i}/alldat.1.${i}_${j}.$r.dat > ./tmp/raw.${r}.dat;
    ../newbar/newbar1 -f ./tmp/raw.${r}.dat | tee -a ./tmp/full.${r}.output

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
awk '/da_est/ {sum+=$2; i++; printf("%3d%10.3f%10.3f\n",i,$2,sum)}' ./tmp/full.${r}.output | tee ./data/full.${r}.dat

