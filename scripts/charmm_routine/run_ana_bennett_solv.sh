idxmax=$1
r=${2:-1}
method=${3:-1}

mkdir -p ./tmp
mkdir -p ./data
rm -f ./tmp/full.${r}.output
#touch tmp/full.result tmp/part1.result tmp/part2.result tmp/part3.result tmp/part4.result #part5.result

for (( i=1; i<${idxmax}; i++ ))
do   
    j=$(($i+1));
    echo "$i -> $j"
    paste intst${i}/prod_${method}.${r}.${i}_${i}.dat intst${j}/prod_${method}.${r}.${j}_${i}.dat \
	  intst${j}/prod_${method}.${r}.${j}_${j}.dat intst${i}/prod_${method}.${r}.${i}_${j}.dat > ./tmp/raw.${r}.dat;
    ../newbar/newbar1 -f ./tmp/raw.${r}.dat | tee -a ./tmp/full.${r}.output

done

echo ""
echo "== Cum. dA =="
echo 
awk '/da_est/ {sum+=$2; i++; printf("%3d%10.3f%10.3f\n",i,$2,sum)}' ./tmp/full.${r}.output | tee ./data/full.${r}.dat

