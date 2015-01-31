SCRIPTS=/home/zchads1/scripts
PEDIGREE=pedigree.pro
python $SCRIPTS/generate_hapdrop.py \
SCRIPTS=/home/zchads1/scripts
PEDIGREE=pedigree.pro
python $SCRIPTS/generate_hapdrop.py \
$PEDIGREE \
cm_lengths.txt \
0.1 31 8


GENOTYPES=/home/zchads1/cluster/exome/annotate/AJ/AJ_moderate.gen
MERLIN=~/merlin/merlin
OUT=.

i=1

$MERLIN -p for_simulation.ped \
-d for_simulation.dat \
-m for_simulation.map \
-f for_simulation.frq \
--bits 0 --simulate -r ${i} \
--save --prefix $OUT/sim_${i} \
--swap --megabytes 1000 > $OUT/run_${i}.log

python $SCRIPTS/extract_from_simulated_multiple_v2.py \
$OUT/sim_${i}-replicate.ped \
$OUT/sim_${i}-replicate.dat \
8 \
$OUT/sim_flow_${i} >> $OUT/run_${i}.log

python ../hapdrop_v5_unphased.py \
pedigree.pro \
../AJ_moderate.cm \
$OUT/sim_flow_${i}_maternal.txt \
$OUT/sim_flow_${i}_paternal.txt \
$GENOTYPES \
$OUT/${i} >> $OUT/run_${i}.log


