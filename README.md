# Fine_Mapping

Identification of causal variants within genes that have been previously associated with Parkinson's Disease.\
All data and results can viewed and downloaded via the following interactive Rmarkdown output files below:  

### Current Results
* [Fine mapping with Parkinson's and Alzheimer's Disease GWAS Summary Stat Data](https://rajlabmssm.github.io/Fine_Mapping/Fine_Mapping.html)


### Alternative Analyses
* [Using LD matrix calculated from 1000 Genomes Project Phase 1 data (1/7/13); European superpopulation only](https://rajlabmssm.github.io/Fine_Mapping/Fine_Mapping_1KGphase1_EUR.html)
* [Using LD matrix calculated from 1000 Genomes Project Phase 3 data (7/14/16); European superpopulation only](https://rajlabmssm.github.io/Fine_Mapping/Fine_Mapping_1KGphase3_EUR.html)
* [Using LD matrix calculated from 1000 Genomes Project Phase 3 data (7/14/16); All populations](https://rajlabmssm.github.io/Fine_Mapping/Fine_Mapping_1KGphase3_allPops.html)


## Methods

We utilize a variety of statistical and functional fine mapping methodologies, including:
* [susieR](https://github.com/stephenslab/susieR)
* [CAVIAR](http://genetics.cs.ucla.edu/caviar/)
* [DAP](https://github.com/xqwen/dap)

## Data

* The PD GWAS summary statistics used here are from the [Nalls et al. (2018) preprint](https://www.biorxiv.org/content/10.1101/388165v1). Summary statistics (limited to only the top SNPs identified by Nalls et al. (2018) via a different fine mapping protocol) can be [found here](https://github.com/neurogenetics/meta5)
* LD matrices are calculated using .vcf files from the [1000 Genomes Project](http://www.internationalgenome.org/) (Phase 1 or Phase 3), and the LD function within the R package [gaston](https://cran.r-project.org/web/packages/gaston/gaston.pdf)
* Both Nalls et al. (2019) summary stats and 1000 Genomes Project LD calculations used human genome annotation GRCh37


## Author

<a href="https://bschilder.github.io/BMSchilder/" target="_blank">Brian M. Schilder, Bioinformatician II</a>\
<a href="https://rajlab.org" target="_blank">Raj Lab</a>\
<a href="https://icahn.mssm.edu/about/departments/neuroscience" target="_blank">Department of Neuroscience, Icahn School of Medicine at Mount Sinai</a>\

