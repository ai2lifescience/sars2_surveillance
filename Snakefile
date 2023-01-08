import os,sys,time,shutil,json
configfile: 'sample_config/test.yaml'


##############################
######## Init
##############################
Lsample=config['Lsample']
batch_name=config['batch_name']

# structure
indir=config['indir']
outdir=config['outdir']
Plog=f'{outdir}/0log'
Pstat=f'{outdir}/0stat'
Pqc=f'{outdir}/1qc'
Pmap=f'{outdir}/2map'
Pdedup=f'{outdir}/3dedup'
Pconsensus=f'{outdir}/4consensus'
Plineage=f'{outdir}/5lineage'
Ptree=f'{outdir}/6tree'


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
        smp_genome_fa = expand(Pmap + '/{sample}/{sample}_smp_genome.fa', sample=Lsample),
        sort_bam = expand(Pmap + '/{sample}/{sample}_sort.bam', sample=Lsample),
        rmp_sort_bam = expand(Pdedup + '/{sample}/{sample}_rmPrimer_sort.bam', sample=Lsample),
        dedup_bam = expand(Pdedup + '/{sample}/{sample}_dedup.bam', sample=Lsample),
        mask_genome_fa = expand(Pconsensus + '/{sample}/{sample}_mask_genome.fa', sample=Lsample),
        consensus_fa = expand(Pconsensus + '/{sample}/{sample}_consensus.fa', sample=Lsample),
        lineage_csv = expand(Plineage + '/{sample}/{sample}_lineage.csv', sample=Lsample),
        upstream_stat = expand(Pstat + '/{batch_name}_upstream_stat.csv', batch_name=batch_name),
        merge_consensus_fa = expand(Pstat + '/{batch_name}_consensus.fa', batch_name=batch_name),
        nextclade_info = expand(Ptree + '/{batch_name}_nextclade.csv',batch_name=batch_name),
        tree_json = expand(Ptree + '/{batch_name}.auspice.json',batch_name=batch_name),
        index = expand(Ptree + '/{batch_name}_sequence_index.tsv',batch_name=batch_name),
        alignment = expand(Ptree + '/{batch_name}_aligned.fasta',batch_name=batch_name),
        tree = expand(Ptree + '/{batch_name}_tree_raw.nwk',batch_name=batch_name)   


##################################
### Upstream
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
    resources: cpus=config['qc_cpus']
    params: 
        fq_cut_mean_qual=config['fq_cut_mean_qual'],
        fq_min_len=config['fq_min_len']
    conda: 'envs/surveillance.yaml'
    shell:"""
        fastp --thread {resources.cpus} -i {input.r1} -I {input.r2} -h {output.html} -j {output.json} \\
            -o {output.trim_r1} -O {output.trim_r2} --unpaired1 {output.unpair_r1} --unpaired2 {output.unpair_r2} \\
            -q 15 -u 40 -l {params.fq_min_len} --cut_right --cut_window_size 20 --cut_mean_quality {params.fq_cut_mean_qual} --correction \\
            1>{log.o} 2>{log.e}
        """

rule map:
    input: 
        trim_r1 = rules.qc.output.trim_r1,
        trim_r2 = rules.qc.output.trim_r2,
        genome_fa = config['genome_fa']
    output:
        smp_genome_fa = Pmap + '/{sample}/{sample}_smp_genome.fa',
        sam = Pmap + '/{sample}/{sample}.sam',
        bam = Pmap + '/{sample}/{sample}.bam',
        sort_bam = Pmap + '/{sample}/{sample}_sort.bam',
    log: e = Plog + '/map/{sample}.e', o = Plog + '/map/{sample}.o'
    benchmark: Plog + '/map/{sample}.bmk'
    resources: cpus=config['map_cpus']
    conda: 'envs/surveillance.yaml'
    shell:"""
        # create sample specific genome fasta
        grep ">" {input.genome_fa} | sed "s/>.*/>{wildcards.sample}/g" > {output.smp_genome_fa}
        awk 'NR>1{{printf "%s",$0}}' {input.genome_fa} >> {output.smp_genome_fa}
        
        # map by bwa
        # bwa index {output.smp_genome_fa} 1>{log.o} 2>{log.e}
        # bwa mem -t {resources.cpus} -v 1 {output.smp_genome_fa} \\
        # {input.trim_r1} {input.trim_r2} -o {output.sam} 1>>{log.o} 2>>{log.e}
        # sambamba view --sam-input -o {output.bam} -f bam -t {resources.cpus} {output.sam} 1>>{log.o} 2>>{log.e}
        # sambamba sort -n -o {output.sort_bam} -t {resources.cpus} {output.bam} 1>>{log.o} 2>>{log.e}
        
        # map by minimap2
        minimap2 -t {resources.cpus} -ax sr {output.smp_genome_fa} {input.trim_r1} {input.trim_r2} -o {output.sam} 1>{log.o} 2>{log.e}
        sambamba view --sam-input -o {output.bam} -f bam -t {resources.cpus} {output.sam} 1>>{log.o} 2>>{log.e}
        sambamba sort -n -o {output.sort_bam} -t {resources.cpus} {output.bam} 1>>{log.o} 2>>{log.e}
        
        """

rule rm_primer:
    input: 
        sort_bam = rules.map.output.sort_bam,
        primer_bed = config['primer_bed']
    output:
        fix_bam = Pdedup + '/{sample}/{sample}_fixmate.bam',
        fix_sort_bam = Pdedup + '/{sample}/{sample}_fixmate_sort.bam',
        rmp_sort_bam = Pdedup + '/{sample}/{sample}_rmPrimer_sort.bam',
    log: e = Plog + '/rm_primer/{sample}.e', o = Plog + '/rm_primer/{sample}.o'
    benchmark: Plog + '/rm_primer/{sample}.bmk'
    resources: cpus=config['rm_primer_cpus']
    params: 
        tmp=Pdedup + '/{sample}/{sample}_tmp',
        rmp_bam_prefix=Pdedup + '/{sample}/{sample}_rmPrimer'
    conda: 'envs/surveillance.yaml'
    shell:"""
    	samtools fixmate -m {input.sort_bam} {output.fix_bam} -@ {resources.cpus} 1>>{log.o} 2>>{log.e}
        sambamba sort -o {output.fix_sort_bam} --tmpdir {params.tmp} -t {resources.cpus} {output.fix_bam} 1>>{log.o} 2>>{log.e}
        ivar trim -i {output.fix_sort_bam} -b {input.primer_bed} -p {params.rmp_bam_prefix} -e 1>>{log.o} 2>>{log.e}
        sambamba sort -o {output.rmp_sort_bam} --tmpdir {params.tmp} -t {resources.cpus} {params.rmp_bam_prefix}.bam 1>>{log.o} 2>>{log.e}
        """

rule dedup:
    input: rmp_sort_bam = rules.rm_primer.output.rmp_sort_bam
    output: dedup_bam = Pdedup + '/{sample}/{sample}_dedup.bam'
    log: e = Plog + '/dedup/{sample}.e', o = Plog + '/dedup/{sample}.o'
    benchmark: Plog + '/dedup/{sample}.bmk'
    resources: cpus=config['dedup_cpus']
    conda: 'envs/surveillance.yaml'
    shell:"""
        sambamba markdup -r -t {resources.cpus} --overflow-list-size=500000 {input.rmp_sort_bam} {output.dedup_bam} 1>>{log.o} 2>>{log.e}
        sambamba index -t {resources.cpus} {output.dedup_bam} 1>>{log.o} 2>>{log.e}
        """

rule mask_genome:
    input: 
        dedup_bam = rules.dedup.output.dedup_bam,
        smp_genome_fa = rules.map.output.smp_genome_fa
    output: 
        bedgraph = Pconsensus + '/{sample}/{sample}.bedgraph',
        lowcov_bed = Pconsensus + '/{sample}/{sample}_lowcov.bed',
        mask_genome_fa = Pconsensus + '/{sample}/{sample}_mask_genome.fa'
    log: e = Plog + '/mask_genome/{sample}.e', o = Plog + '/mask_genome/{sample}.o'
    benchmark: Plog + '/mask_genome/{sample}.bmk'
    resources: cpus=config['mask_genome_cpus']
    params: min_coverage=config['min_coverage']
    conda: 'envs/surveillance.yaml'
    shell:"""
        bedtools genomecov -ibam {input.dedup_bam} -bga > {output.bedgraph} 2>>{log.e}
        awk -v cov={params.min_coverage} '$4<cov' {output.bedgraph} | bedtools merge -i - >{output.lowcov_bed} 2>>{log.e}
        bedtools maskfasta -fi {input.smp_genome_fa} -bed {output.lowcov_bed} -fo {output.mask_genome_fa} 1>>{log.o} 2>>{log.e}
        """

rule consensus:
    input: 
        dedup_bam = rules.dedup.output.dedup_bam,
        smp_genome_fa = rules.map.output.smp_genome_fa,
        mask_genome_fa = rules.mask_genome.output.mask_genome_fa
    output: 
        bqsr_bam = Pconsensus + '/{sample}/{sample}_bqsr.bam',
        pileup = Pconsensus + '/{sample}/{sample}.pileup',
        variant_info = Pconsensus + '/{sample}/{sample}_variant_info.tsv',
        consensus_fa = Pconsensus + '/{sample}/{sample}_consensus.fa'
    log: e = Plog + '/consensus/{sample}.e', o = Plog + '/consensus/{sample}.o'
    benchmark: Plog + '/consensus/{sample}.bmk'
    resources: cpus=config['consensus_cpus']
    params: 
        consensus_prefix=Pconsensus + '/{sample}/{sample}_consensus',
        variant_info_prefix=Pconsensus + '/{sample}/{sample}_variant_info',
        min_coverage=config['min_coverage'],
        min_allele_freq=config['min_allele_freq']
    conda: 'envs/surveillance.yaml'
    shell:"""
        lofreq indelqual {input.dedup_bam} --dindel -f {input.mask_genome_fa} -o {output.bqsr_bam} 1>>{log.o} 2>>{log.e}
        samtools index {output.bqsr_bam} -@ {resources.cpus} 1>>{log.o} 2>>{log.e}
        # mpileup
        samtools mpileup -aa -A -d 0 -Q 0 {output.bqsr_bam} \\
            -o {output.pileup} 1>>{log.o} 2>>{log.e}
        # variant info
        cat {output.pileup} | ivar variants -p {params.variant_info_prefix} \\
            -r {input.smp_genome_fa} \\
            -t {params.min_allele_freq} -m {params.min_coverage} 1>>{log.o} 2>>{log.e}
        # consensus
        cat {output.pileup} | ivar consensus -p {params.consensus_prefix} \\
            -t {params.min_allele_freq} -m {params.min_coverage} 1>>{log.o} 2>>{log.e}
        sed -i "s/^>.*/>{wildcards.sample}/g" {output.consensus_fa} 1>>{log.o} 2>>{log.e}
        """

##################################
### Downstream
##################################

rule pangolin:
    input: consensus_fa = rules.consensus.output.consensus_fa
    output: lineage_csv = Plineage + '/{sample}/{sample}_lineage.csv',
    log: e = Plog + '/pangolin/{sample}.e', o = Plog + '/pangolin/{sample}.o'
    benchmark: Plog + '/pangolin/{sample}.bmk'
    resources: cpus=config['pangolin_cpus']
    params: pangolin_prifix=Plineage + '/{sample}/{sample}',
    conda: 'envs/surveillance.yaml'
    shell:"""
        pangolin {input.consensus_fa} \\
            --outfile {output.lineage_csv} -t {resources.cpus} 1>>{log.o} 2>>{log.e}
        """

rule notebook_init:
    input: 
        upstream_stat_r='scripts/upstream_stat.r.ipynb'
    output: 
        upstream_stat_r=Plog + '/upstream_stat.r.ipynb'
    resources: cpus=1
    log: e = Plog + '/notebook_init.e', o = Plog + '/notebook_init.o'
    run:
        for k, f in input.items():
            notebook = json.loads(open(f, 'r').read())
            if f.endswith('.r.ipynb'):
                notebook['metadata']['kernelspec']['name'] = 'ir'
                notebook['metadata']['kernelspec']['display_name'] = 'R'
            elif f.endswith('.py.ipynb'):
                notebook['metadata']['kernelspec']['name'] = 'python3'
                notebook['metadata']['kernelspec']['display_name'] = 'Python 3 (ipykernel)'
            with open(output[k], 'w') as out:
                json.dump(notebook, out)

rule upstream_stat:
    input:
        Lqc = expand(rules.qc.output.json, sample=Lsample),
        Lconsensus = expand(rules.consensus.output.consensus_fa, sample=Lsample),
        Ldup_rate = expand(rules.dedup.log.e, sample=Lsample),
        Lalign_rate = expand(rules.rm_primer.log.o, sample=Lsample),
        Lvariant_info = expand(rules.consensus.output.variant_info, sample=Lsample),
        upstream_stat_r = rules.notebook_init.output.upstream_stat_r
    output: 
        upstream_stat = Pstat + '/{batch_name}_upstream_stat.csv',
        merge_consensus_fa = Pstat + '/{batch_name}_consensus.fa'
    log: 
        notebook = Plog + '/upstream_stat/{batch_name}.r.ipynb', 
        e = Plog + '/upstream_stat/{batch_name}.e', 
        o = Plog + '/upstream_stat/{batch_name}.o'
    benchmark: Plog + '/upstream_stat/{batch_name}.bmk'
    resources: cpus=config['upstream_stat_cpus']
    params: high_altfreq_threshold=config['high_altfreq_threshold']
    conda: 'envs/jupyter.yaml'
    notebook: rules.notebook_init.output.upstream_stat_r


rule nextclade:
    message:
        '''
        Calling mutations
        Adding consensus to the existing tree
        '''
    input: consensus_fa = rules.upstream_stat.output.merge_consensus_fa
    output: 
        nextclade_info = Ptree + '/{batch_name}_nextclade.csv',
        tree_json = Ptree + '/{batch_name}.auspice.json',
    log:  e = Plog + '/nextclade/{batch_name}.e', o = Plog + '/nextclade/{batch_name}.o'
    benchmark: Plog + '/nextclade/{batch_name}.bmk'
    resources: cpus=config['nextclade_cpus']
    params: built_tree=config['built_tree'],
    conda: 'envs/buildTree.yaml'
    shell:
        """
        nextclade run \
            --in-order \
            --input-dataset={params.built_tree} \
            --output-csv={output.nextclade_info} \
            --output-tree={output.tree_json} \
            {input.consensus_fa}
        """

rule augur:
    message:
        '''
        Building a tree with only the consensus 
        '''
    input: consensus_fa = rules.upstream_stat.output.merge_consensus_fa
    output: 
        index = Ptree + '/{batch_name}_sequence_index.tsv',
        alignment = Ptree + '/{batch_name}_aligned.fasta',
        tree = Ptree + '/{batch_name}_tree_raw.nwk',
    log:  e = Plog + '/augur/{batch_name}.e', o = Plog + '/augur/{batch_name}.o'
    benchmark: Plog + '/augur/{batch_name}.bmk'
    resources: cpus=config['augur_cpus']
    params: 
        reference = config['root'],
    conda: 'envs/buildTree.yaml'
    shell:
        """
        augur index \
            --sequences {input.consensus_fa} \
            --output {output.index}
        augur align \
            --sequences {input.consensus_fa} \
            --reference-sequence {params.reference} \
            --output {output.alignment} \
            --fill-gaps
        augur tree \
            --alignment {output.alignment} \
            --nthreads {resources.cpus} \
            --output {output.tree} 
        """

