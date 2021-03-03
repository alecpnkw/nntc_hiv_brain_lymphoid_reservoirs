# LP71-199_hiv_brain_lymphoid_reservoirs

### Overview

Sequence post-processing of NNTC env sequences from brain and lymphoid tissues. SGA consensus sequences used as input to the workflow are found at `resources/sequences.fa`. Output from the geno2pheno coreceptor [tool](https://coreceptor.geno2pheno.org/) can be found at `config/metadata`. 

This workflow performs functional filtering of env sequences, multiple alignment with mafft, and hypermutation detection via a python [implementation](https://gist.github.com/philiptzou/8d6c7c61d2242f730a1a2f87ba9a2a72) by @philiptzou of the LANL tool [Hypermut2](https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html). 

FST and diversity analysis by participant and tissue are contained in this [notebook](workflow/notebooks/Fst_diversity_stat.ipynb)


### Setup and Usage

Install [snakemake](https://snakemake.readthedocs.io/en/stable/) and [biopython](https://biopython.org/wiki/Packages) using conda. 

```
conda install -c conda-forge snakemake-minimal
conda install -c conda-forge biopython
```

Install Julia > 1.2 from [here](https://julialang.org/downloads/) and initialize the Julia environment in the project directory. Next, run

```
snakemake -n --use-conda
```

to preview jobs and

```
snakemake --jobs {n} --use-conda
```

to execute, where "n" is the number of parallel jobs desired. 
