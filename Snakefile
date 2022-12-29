import os,sys,time,shutil
configfile: 'sample_config/test.yaml'


##############################
######## Init
##############################
Lsample=config['Lsample']
Dresources=config['Dresources']

# ref
# ref_genome = config['Dref_genome'][genome]

# structure
indir=config['indir']
outdir=config['outdir']
Plog=f'{outdir}/0log'
Pqc=f'{outdir}/1qc'
Pmap=f'{outdir}/2map'



for p in [Pqc,Pmap]: 
    os.makedirs(p, exist_ok=True)

onstart:
    for r in workflow.rules:
        os.makedirs(f'{Plog}/{r}', exist_ok=True)
    global time_start, smk_file_name
    time_start = time.strftime("%y%m%d_%H%M%S", time.localtime())
    smk_file_name = workflow.snakefile.split('/')[-1]
    shutil.copyfile(workflow.snakefile, f'{Plog}/all/{time_start}_{smk_file_name}')

##################################
### Rule all
##################################
rule all:
    input:
        trim_r1 = expand(Pqc + '/{sample}/{sample}_trim_R1.fastq.gz', sample=Lsample),
        trim_r2 = expand(Pqc + '/{sample}/{sample}_trim_R2.fastq.gz', sample=Lsample),


##################################
### Common rules
##################################
rule qc:
    input: 
        r1 = indir + '/{sample}_R1.fastq.gz',
        r2 = indir + '/{sample}_R2.fastq.gz'
    output:
        trim_r1 = Pqc + '/{sample}/{sample}_trim_R1.fastq.gz',
        trim_r2 = Pqc + '/{sample}/{sample}_trim_R2.fastq.gz',
        unpair_r1 = Pqc + '/{sample}/{sample}_unpair_R1.fastq.gz',
        unpair_r2 = Pqc + '/{sample}/{sample}_unpair_R2.fastq.gz',
        html = Pqc + '/{sample}/{sample}.html',
        json = Pqc + '/{sample}/{sample}.json',
    log: e = Plog + '/qc/{sample}.e', o = Plog + '/qc/{sample}.o'
    benchmark: Plog + '/qc/{sample}.bmk'
    resources: cpus=Dresources['qc_cpus']
    conda: 'envs/surveillance.yml'
    shell:"""
        fastp --thread {resources.cpus} -i {input.r1} -I {input.r2} -h {output.html} -j {output.json} \\
            -o {output.trim_r1} -O {output.trim_r2} --unpaired1 {output.unpair_r1} --unpaired2 {output.unpair_r2} \\
            -q 15 -u 40 -l 25 --cut_right --cut_window_size 20 --cut_mean_quality 30 --correction \\
            1>{log.o} 2>{log.e}
        """
