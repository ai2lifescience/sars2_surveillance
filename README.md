# sars2_surveillance
Variants and Genomic Surveillance for SARS-CoV-2


## run pipeline test via snakemake

1. install snakemake
```sh
conda install snakemake
```
2. modified the `sample_config/test.yaml` if you want
3. go back to the root directory where `Snakefile` in
```sh
snakemake -j 5
```
- `j`: the number of jobs to run in parallel
- the first run will create an anonymous conda environment, which can be found by `conda env list`
- the output of test data will be in `test/output`
