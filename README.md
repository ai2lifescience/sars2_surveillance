# sars2_surveillance
Variants and Genomic Surveillance for SARS-CoV-2


## run pipeline test via snakemake

1. Install snakemake
```sh
conda install -c conda-forge -c bioconda snakemake
```
or
```sh
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```
Strict channel priority can dramatically speed up conda operations and also reduce package incompatibility problems.
```sh
conda config --set channel_priority strict
```
2. Modify the `sample_config/test.yaml` to customize the pipeline parameters if you want
3. Go back to the root directory where `Snakefile` in
```sh
snakemake -j 5 --use-conda
```
or
```sh
snakemake -j 5 --use-conda --configfile your_configfile_path/new_batch.yaml
```
- `j`: the number of jobs to run in parallel
- The first run will create two anonymous conda environments, one of which includes softwares in `envs/surveillance.yaml`, the other includes softwares in `envs/jupyter.yaml`. So you don't need any prepared environment
- The anonymous conda environments can be found by `conda env list`
- The output of test data will be in `test/output`
4. If you want to run on a computation server
- For `slurm` system, to assign 4g memory for each CPU and run on `compute` partition, use
```sh
snakemake -j 5 --use-conda --cluster "sbatch -p compute --mem-per-cpu 4g -c {resources.cpus} -J {rule}_{wildcards} -o {log.o} -e {log.e}"
```
- This command will generate a slurm job for each rule-sample combination
- For convenience, you can write a config file in `~/.config/snakemake`, for example
```sh
mkdir -p ~/.config/snakemake/my_smk_config
echo 'cluster: "sbatch -p compute --mem-per-cpu 4g -c {resources.cpus} -J {rule}_{wildcards} -o {log.o} -e {log.e}"' > ~/.config/snakemake/my_smk_config/config.yaml
echo 'use-conda: TRUE' >> ~/.config/snakemake/my_smk_config/config.yaml
```
- Then you can run pipeline by
```sh
snakemake -j 5 --profile my_smk_config
```
5. Update nextclade dataset
```sh
 nextclade dataset get --name 'sars-cov-2' --output-dir 'src/sars-cov-2'
```

## run pipeline for new sample batch

1. Rename the rawdata into `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz` just like the tests
2. Copy the contents of `sample_config/test.yaml`, modify `indir`, `outdir`, `Lsample` (the sample names of this run batch), `batch_name` and so on. Save it into another file, `sample_config/new_batch.yaml` for example
3. Assign the config file and run the pipeline by
```sh
snakemake -j 5 --profile my_smk_config --configfile sample_config/new_batch.yaml
```

## Credit and Acknowledgements
- The script to convert variants.tsv files into .vcf files was obtained from: https://github.com/nf-core/viralrecon/blob/dev/bin/ivar_variants_to_vcf.py