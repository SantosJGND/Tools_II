## Imputation repository.


This sub-directory holds work on imputation. As usual, this research is done with population genetic data in mind. 

The data consists a single genotype data set. Variables are variant count features ranging between 0 and 2; Samples are designed to derive from a semi-consistent population network. _Semi-consistent_ is used here to indicate that certain observations have variable pdfs, and the characteristics of the structure vary (cluster distance may change).

### Data generation

VCF files are generated using the [Genome Simulator](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Simulate_genomes/Genomic%20structure%20Simulator.ipynb) tool of the first Tools repository [link](https://github.com/SantosJGND/Tools_and_toys).

- replicated here for the specific data sets used [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/prepare_vcfs.ipynb).


### I. Distances / Dimensionality reduction. 

Window based analysis constructs data sets of distance data with which to predict position of missing observation in incomplete data set.

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_II/blob/master/Imputation/Impute_I_distances.ipynb)

