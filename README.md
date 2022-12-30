# sars2_surveillance
Variants and Genomic Surveillance for SARS-CoV-2


## run pipeline test via snakemake

1. Install snakemake
```sh
conda install snakemake
```
2. Modify the `sample_config/test.yaml` to customize the pipeline parameters if you want
3. Go back to the root directory where `Snakefile` in
```sh
snakemake -j 5 --use-conda
```
- `j`: the number of jobs to run in parallel
- The first run will create an anonymous conda environment, which includes all the software needed according to `envs/surveillance.yaml`. So you don't need any prepared environment.
- The anonymous conda environment can be found by `conda env list`
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
