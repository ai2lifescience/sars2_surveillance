# sars2_surveillance
Variants and Genomic Surveillance for SARS-CoV-2


## run pipeline test via snakemake

1. Install snakemake
```sh
conda install snakemake
```
2. Modify the `sample_config/test.yaml` to customize the pipeline parameters if you want
3. Go back the root directory where `Snakefile` in
```sh
snakemake -j 5
```
- `j`: the number of jobs to run in parallel
- The first run will create an anonymous conda environment, which includes all the software needed according to `envs/surveillance.yaml`. So you don't need any prepared environment.
- The anonymous conda environment can be found by `conda env list`
- The output of test data will be in `test/output`
