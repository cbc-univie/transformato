# concatenate indiv. sets and analyse full data set

states=$1 # number intst dirs - 1
#r=${2:-1}
#blk=10000

mkdir -p ./tmp
rm -f ./tmp/gasp.cum.output
rm -f ./tmp/raw*dat

for (( i=1; i<=${states}; i++ ))
do
    j=$(($i+1));
    echo "$i -> $j"
    cat intst${i}/gasp.?.${i}_${i}.dat > intst${i}/gasp.cum.${i}_${i}.dat
    cat intst${i}/gasp.?.${i}_${j}.dat > intst${i}/gasp.cum.${i}_${j}.dat
    cat intst${j}/gasp.?.${j}_${i}.dat > intst${j}/gasp.cum.${j}_${i}.dat
    cat intst${j}/gasp.?.${j}_${j}.dat > intst${j}/gasp.cum.${j}_${j}.dat
    paste intst${i}/gasp.cum.${i}_${i}.dat intst${j}/gasp.cum.${j}_${i}.dat \
	  intst${j}/gasp.cum.${j}_${j}.dat intst${i}/gasp.cum.${i}_${j}.dat > ./tmp/raw.cum.$i.dat;
    ~/work/festutorial/absolute/newbar/newbar1 -f ./tmp/raw.cum.$i.dat | tee -a ./tmp/gasp.cum.output

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
awk '/da_est/ {sum+=$2; i++; printf("%3d%10.3f%10.3f\n",i,$2,sum)}' ./tmp/gasp.cum.output | tee ./data/gasp.cum.dat

