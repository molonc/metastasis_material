
# Metastasis manuscript materials

This repository contains all materials and scripts used in manuscript: 
```
Hoa Tran, Gurdeep Singh, Hakwoo Lee, Damian Yap, Eric Lee, William Daniels, Farhia Kabeer, Ciara H O’Flanagan, Vinci Au, Michael Van Vliet, Daniel Lai, Elena Zaikova, Sean Beatty, Esther Kong, Shuyu Fan, Jessica Chan,Sam Dang, Viviana Cerda, Teresa Ruiz de Algaza, Andrew Roth, Samuel Aparicio

Somatic copy number mutations contribute to fitness in transplantation models of spontaneous human breast cancer metastasis.
biorxiv preprint DOI: https://www.biorxiv.org/content/ (will update link soon)

```

- [Overview](#overview)
- [Materials](#materials)
  - [Uploaded data link](#uploaded-data-link)
  - [Meta data](#meta-data)
  - [Phylogenetic tree results](#phylogenetic-tree-results)
  - [Tissue screening](#tissue-screening)
  - [Main figures](#main-figures)
  - [Supplementary figures](#supplementary-figures)
  - [Supplementary tables](#supplementary-tables)
 
- [Scripts](#scripts)
  - [Copy number analysis](#copy-number-analysis)
  - [Preprocessing functions](#preprocessing-functions)
  - [Differential expression analysis](#differential-expression-analysis)
  
- [Citation](#citation)
- [Licence](#licence)


## Overview

The contribution of somatic gene dosage mutations (CNA) to breast cancer metastasis remains poorly defined. Using 9 transplantable human triple-negative breast cancer xenografts, we studied the fitness of copy number clones in spontaneous metastasis from orthotopic transplant sites. Metastatic site preference was strongly patient-dependent, and the emergence of metastases exhibited a general trend toward slower growth at the orthotopic site. In our models, single-cell whole-genome sequencing of primary and metastatic sites showed that distant metastases were most often the result of minor prevalence clones at the orthotopic site. This suggests some metastatic phenotypes may be weakly negatively fit at the primary site. We validated the existence of a fitness hierarchy of copy number clones using a previously established paradigm of re-mixing and re-transplanting clones. Single-cell clone analysis of competitive repopulation and re-emergence of metastases showed that CNAs arising in cancer evolution can mediate metastatic fitness. Moreover, some clones exhibiting strong metastatic tendency exhibited weaker survival at the primary site compared with clones of lesser metastatic potential. This is consistent with the notion that metastatic phenotypes could have a fitness cost at the primary site. Finally, we conducted RNA-seq analysis combined with DriverNet analysis to dissect the contribution of CNA-mediated versus genome-independent transcriptional states. CNA mutations appeared to contribute strongly to transcriptional differences between clones. Among clones of high metastatic potential, we observed CNA-mediated and CNA-independent convergence on pathways such as EMT, established as mediators of metastatic cell survival at distant sites. Taken together, our data point to a contribution of CNA-mediated cancer evolution to the metastatic states and identify distant site context as a key determinant of CNA-mediated fitness. 




## Materials

### Uploaded data link
- DLP+ alignment copy number data libraries using HMMCopy method are at: [Uploaded Data URL](https://ega-archive.org/studies/EGAS00000000000) - in editing
- bulk RNA-seq Kallisto alignment libraries are at: [Uploaded Data URL](https://ega-archive.org/studies/EGAS00000000000) - in editing


### Meta data
- Library infos, primary and metastasis statuses and time series - passages are noted at [materials/dlp_trees](https://github.com/molonc/metastasis_material/tree/main/materials/dlp_trees/EGA_metastasis_project/).
- The full details of metadata files are noted at [materials/inventory_19Mar2024](https://github.com/molonc/metastasis_material/tree/main/materials/inventory_19Mar2024)


### Phylogenetic tree results

- The inferred Sitka phylogeny trees for patient Pt1-SA919, Pt2-SA535 [materials/dlp_trees/](https://github.com/molonc/metastasis_material/tree/main/materials/dlp_trees/)
- The inferred Sitka phylogeny trees for patient Pt1-SA919 at mixture experiments [materials/dlp_trees/SA919_mixing_experiment/](https://github.com/molonc/metastasis_material/tree/main/materials/dlp_trees/SA919_mixing_experiment)



### Tissue screening
- TMA score results for TMA microarray tissue screening are at file  [materials/tumour_volumes/](https://github.com/molonc/metastasis_material/tree/main/materials/tumour_volumes)

### Cis trans genes 
- DE analysis are at [materials/bulkRNAseq/](https://github.com/molonc/metastasis_material/tree/main/materials/bulkRNAseq/) 
- Pathway analysis results are at [materials/bulkRNAseq/](https://github.com/molonc/metastasis_material/tree/main/materials/bulkRNAseq/)

### Main figures 
- Figure files at [materials/main_figures/](https://github.com/molonc/metastasis_material/tree/main/main_figures/) 


### Supplementary figures 
- Supplementary figure files at [materials/supplementary_figures/](https://github.com/molonc/metastasis_material/tree/main/supplementary_figures/) 

### Supplementary tables
- Supplementary tables at [materials/supplementary_tables/](https://github.com/molonc/metastasis_material/tree/main/supplementary_tables/) 


## Scripts
### Copy number analysis
- Script for phylogenetic tree reconstruction, clonal analyses are at [scripts/corrupt_tree/src](https://github.com/molonc/metastasis_material/tree/main/scripts/corrupt_tree/src/)

### Preprocessing functions
- Preprocessing scripts for bulk RNA-seq analysis are at [scripts/bulk_rna/](https://github.com/molonc/metastasis_material/tree/main/scripts/bulk_rna/)

### Differential expression analysis
- Scripts are at [scripts/bulk_rna/](https://github.com/molonc/metastasis_material/tree/main/scripts/bulk_rna/)




## Citation
```
Hoa Tran, Gurdeep Singh, Hakwoo Lee, Damian Yap, Eric Lee, William Daniels, Farhia Kabeer, Ciara H O’Flanagan, Vinci Au, Michael Van Vliet, Daniel Lai, Elena Zaikova, Sean Beatty, Esther Kong, Shuyu Fan, Jessica Chan,Sam Dang, Viviana Cerda, Teresa Ruiz de Algaza, Andrew Roth, Samuel Aparicio

Somatic copy number mutations contribute to fitness in transplantation models of spontaneous human breast cancer metastasis.
biorxiv preprint DOI: https://www.biorxiv.org/content/ (will update link soon)

```

## Licence
[Apache licence v2.0](https://github.com/molonc/metastasis_material/blob/main/LICENCE)

