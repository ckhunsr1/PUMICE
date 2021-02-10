<!-- ABOUT THE PROJECT -->
## PUMICE

PUMICE (**P**rediction **U**sing **M**odels **I**nformed by **C**hromatin conformations and **E**pigenomics) is a tool to create gene expression prediction models for transcriptome-wide association studies. Specifically, PUMICE leverages tissue-specific 3D genomic and epigenomic data to define regions that harbor cis-regulatory variants and prioritize them accordingly.

<p align="center">
   <img src="https://github.com/ckhunsr1/PUMICE/blob/master/image/PUMICE_overview.png" width="934" height="166.7">
</p>

<!-- GETTING STARTED -->
## Getting Started

PUMICE requires R 4.0, several R packages, and bedtools.

### Prerequisites

A list of R packages required for PUMICE includes optparse, data.table, tidyr, tidyverse, dplyr, IRanges, GenomicRanges, genefilter, glmnet, caret.


### Tool overview

To run PUMICE, two steps are required.
1. First, we need to run nested cross-validation to determine which window type and penalty factor are optimal (i.e. least mean cross-validated error) for each gene. This step is computationally intensive; therefore, we require users to run this step using parallel computation for the 22 autosomes and each window type. Users can further split each job into multiple jobs using the options total_file_num and file_num. **PUMICE.nested_cv.R** script can be found [here](https://github.com/ckhunsr1/PUMICE/blob/master/Model_training/PUMICE.nested_cv.R).
```
   Rscript PUMICE.nested_cv.R
      --geno [Path to genotype data]
      --chr [Chromosome number]
      --exp [Path to expression data]
      --out [Path to output directory]
      --method [Window type to be used for creating models]
      --type [Specific 3D genome windows being used/Specific constant window size being used (in kb)]
      --window_path [Path to 3D genome window file]
      --bedtools_path [Path to bedtools software]
      --epi_path [Path to epigenomic data]
      --fold [Number of folds to be performed for nested cross-validation]
      --total_file_num [Number of total jobs to be splitted into]
      --file_num [Job number]
      --noclean [Do not delete any temporary files]
   ```
2. Second, we need to run cross-validation to create gene expression prediction model using window type and penalty factor derived from the first step. **PUMICE.compute_weights.R** script can be found [here](https://github.com/ckhunsr1/PUMICE/blob/master/Model_training/PUMICE.compute_weights.R).
```
   Rscript PUMICE.compute_weights.R
      --geno [Path to genotype data]
      --chr [Chromosome number]
      --exp [Path to expression data]
      --out [Path to output directory]
      --pchic_path [Path to pchic window file]
      --loop_path [Path to loop window file]
      --tad_path [Path to tad window file]
      --domain_path [Path to domain window file]
      --bedtool_path [Path to bedtools software]
      --epi_path [Path to epigenomic data]
      --fold [Number of folds to be performed for cross-validation]
      --noclean [Do not delete any temporary files]
   ```

<!-- USAGE EXAMPLES -->
## Usage

We provided example input data [here](https://github.com/ckhunsr1/PUMICE/blob/master/examples/example_input.zip).

Example of shell script used to run step1  can be found [here](https://github.com/ckhunsr1/PUMICE/blob/master/examples/running_nested_cv.sh).
* Outputs from the step1 using example input data are provided [here](https://github.com/ckhunsr1/PUMICE/blob/master/examples/example_output_nestedcv.zip).

Example of shell script used to run step2 can be found [here](https://github.com/ckhunsr1/PUMICE/blob/master/examples/running_compute_weights.sh).
* Outputs from the step2 using example input data are provided [here](https://github.com/ckhunsr1/PUMICE/blob/master/examples/example_output_weights.zip).

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Chachrit (Poom) Khunsriraksakul - [@ChachritK](https://twitter.com/ChachritK) - ckhunsr2@gmail.com



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [Dajiang J. Liu](https://dajiangliu.blog/)

