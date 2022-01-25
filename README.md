# H3N2 HA egg-passaging adaptation


### Dependencies ###
* [python](https://www.python.org/) (version 3.9)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [PEAR](https://github.com/tseemann/PEAR)
* [pandas](https://pandas.pydata.org/)
* [biopython](https://github.com/biopython/biopython)

### Input files ###
* All raw reads in fastq format should be placed in fastq/
* 

### Dependencies installation ###
1. Install dependencies by conda:   
```
conda create -n HA_egg -c bioconda -c anaconda -c conda-forge \
  python=3.9 \
  pear \
  biopython \
  snakemake \
```   

2. Activate conda environment:   
``source activate HA_egg``

### Calling mutations from sequencing data ###
1. Using shell script to call sequence merging and variance calling:   
``./code/pipeline.sh``

2. Call mutation frequency at single aa level:   
``python3 code/count_single_mut_freq.py``

3. Organzie data:   
``python3 code/Merge&organize.py``
    - Output file: [./result/HA_EggPassaging_all_singleAAlevel.xlsx](./result/HA_EggPassaging_all_singleAAlevel.xlsx)
