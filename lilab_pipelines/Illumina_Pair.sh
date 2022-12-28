#PBS -q core24
#PBS -l mem=60gb,nodes=1:ppn=12,walltime=48:00:00
#PBS -o [Dir]/[FileName].map.o
#PBS -e [Dir]/[FileName].map.e
#HSCHED -s SARS2+mapping+varscan
#PPN limit 12
#!/bin/bash
startTime=`date +"%Y-%m-%d %H:%M:%S"`
mkdir -p [Dir]
cd [Dir]
export PATH=/software/biosoft/software/samtools-1.9:$PATH:/software/biosoft/software/fastp:/software/biosoft/software/python/anaconda3-python3-2018/bin:/software/biosoft/software/sortmerna-2.1b:/software/biosoft/software/bwa-0.7.17:/software/biosoft/software/minimap2

#fastp质控去接头
fastp -i [FastaFile1] -I [FastaFile2] -o p1.fq -O p2.fq -l 75 -x -w 12 -q 15 -u 40 --correction --cut_window_size 20 --detect_adapter_for_pe --cut_tail --cut_tail_mean_quality 20 > [FileName].qc 2>&1

#kraken提取冠状病毒reads
/pnas/limk_group/shilsh/kraken2/kraken2 --db /pnas/limk_group/shilsh/kraken2/kraken2_database --threads 12 --output kraken.output --report kraken.report --paired p1.fq p2.fq
/pnas/limk_group/shilsh/python_script/filter.py kraken.output --taxid 2499399 --taxid_col 3 |cut -f 2 > 2499399.reads.id
/xtdisk/limk_group/limk/software/seqtk/seqtk subseq p1.fq 2499399.reads.id > 2499399.R1.fq
/xtdisk/limk_group/limk/software/seqtk/seqtk subseq p2.fq 2499399.reads.id > 2499399.R2.fq

#minimap比对参考基因组
minimap2 -t 12 -ax sr /xtdisk/limk_group/limk/ncov/MN908947.fasta 2499399.R1.fq 2499399.R2.fq | samtools view -@ 12 -b -F 4 -F 2048 -F 256 |samtools sort -@ 12 -o [FileName].2499399.bam

#picard去重复
/software/biosoft/software/jdk1.7.0_17/bin/java -jar /software/biosoft/software/picard-tools-1.119/MarkDuplicates.jar I=[FileName].2499399.bam O=[FileName].2499399.dedup.bam REMOVE_DUPLICATES=true METRICS_FILE=[FileName].2499399.dedup.metrics

#trimBam去除reads两端10bp
/software/biosoft/software/bamUtil/bin/bam trimBam [FileName].2499399.dedup.bam [FileName].2499399.trimed.bam 10 -c
samtools sort -@ 12 [FileName].2499399.trimed.bam -o [FileName].2499399.sorted.trimed.bam

#samtools mpileup统计比对情况
samtools flagstat -@ 12 [FileName].2499399.sorted.trimed.bam > [FileName].2499399.mapped.csv &
samtools mpileup -OsBa --reference /xtdisk/limk_group/limk/ncov/MN908947.fasta -d 3000000 -Q 20 -q 20 [FileName].2499399.sorted.trimed.bam > [FileName].2499399.trimed.mpileup

#varscan产生位点比对结果readcounts文件
/software/biosoft/software/jdk1.7.0_17/bin/java -jar /xtdisk/limk_group/limk/VarScan.v2.3.9.jar readcounts [FileName].2499399.trimed.mpileup --min-coverage 0 --output-file [FileName].2499399.varscan.readcounts

#rm *.bam
#rm *.output
#rm *.report
#rm *.fq
#rm *.id
#rm *.mpileup

endTime=`date +"%Y-%m-%d %H:%M:%S"`
st=`date -d "$startTime" +%s`
et=`date -d "$endTime" +%s`
sumTime=$(($et-$st))
echo "Total time: $sumTime seconds."
wait
