<!-- ABOUT THE PROJECT -->
## PUMICE

PUMICE (Prediction Using Models Informed by Chromatin conformations and Epigenomics) is a tool to create gene expression prediction models for transcriptome-wide association studies. Specifically, PUMICE leverages tissue-specific 3D genomic and epigenomic data to define regions that harbor cis-regulatory variants and prioritize them accordingly.

<!-- GETTING STARTED -->
## Getting Started

PUMICE requires R 4.0 and several R packages.

### Prerequisites

A list of R packages required for PUMICE: optparse, data.table, tidyr, tidyverse, dplyr, IRanges, GenomicRanges, genefilter, glmnet, caret.


### Tool overview

To run PUMICE, two steps are required.
1. First, we need to run nested cross-validation to determine which window type and penalty factor are optimal (i.e. least mean cross-validated error) for each gene.
```
   Rscript PUMICE.nested_cv.R
      --geno [Path to genotype data]
      --chr [Chromosome number]
      --exp [Path to expression data]
      --out [Path to output directory]
      --method [Window type to be used for creating models]
      --type [Specific 3D genome windows being used/Specific constant window size being used (in kb)]
      --window_path [Path to 3D genome windows bed file]
      --bedtool_path [Path to bedtools software]
      --epi_path [Path to epigenomic data]
      --fold [umber of folds to be performed for nested cross-validation]
      --total_file_num [Number of total jobs to be splitted into]
      --file_num [Job number]
      --noclean [Do not delete any temporary files]
   ```
2. Second, we need to run cross-validation to create gene expression prediction model using window type and penalty factor derived from the first step.
```sh
   Rscript PUMICE.compute_weights.R
   ```

<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_



<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/othneildrew/Best-README-Template/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Chachrit (Poom) Khunsriraksakul - [@ChachritK](https://twitter.com/ChachritK) - ckhunsr2@gmail.com



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [Othneil Drew](https://github.com/othneildrew)





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
[product-screenshot]: images/screenshot.png
