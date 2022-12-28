partition=compute_new
cpu=56
mem_per_cpu=4g
tmp=temp_submit
IN=/dshare/xielab/analysisdata/Antibody/songwl/nCoV-evo/rawdata

ls ${IN} | while read smp
# for smp in A1097 A1094 A1083
do

cat << EOF > ${tmp}/${smp}_submit.sh
#!/bin/bash
#SBATCH -J evo_${smp}
#SBATCH -p compute_new
#SBATCH -c ${cpu}
#SBATCH --mem-per-cpu ${mem_per_cpu}
#SBATCH -t 10080
#SBATCH -e log/evo_${smp}.e
#SBATCH -o log/evo_${smp}.o

mkdir -p log
bash HAVoC/HAVoC.sh ${IN}/${smp} \
    1>log/evo_${smp}.o 2>log/evo_${smp}.e

EOF
sbatch ${tmp}/${smp}_submit.sh
done