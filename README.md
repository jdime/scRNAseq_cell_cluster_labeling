# scRNAseq_cell_cluster_labeling

Description
================
This repository contains scripts to run and benchmark scRNA-seq cell cluster labeling methods and is a companion to our paper 'Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data' (Diaz-Mejia JJ et al (2019) [https://www.biorxiv.org/content/10.1101/562082v1].


Scripts
================

**Main wrapper scripts**
---
| Script name |  Task(s) |
| ----------------------------------------------------------- |  ------------------------------------------------ |
| `subsamples_gene_classes_and_runs_enrichment_scripts.R`       |  Main wrapper to run and benchmark cell cluster labeling methods |

**Scripts to run cell type labeling methods**
---
| Script name |  Task(s) |
| -------------------------------------------------------------- |  ------------------------------------------------ |
| `obtains_CIBERSORT_for_MatrixColumns.pl`                       | Runs CIBERSORT using gene expression signatures and a matrix with average gene expressions per gene, per cell cluster) |
| `obtains_GSEA_for_MatrixColumns.pl`                            | Runs GSEA using gene expression signatures and a matrix with average gene expressions per cell cluster) |
| `obtains_GSVA_for_MatrixColumns.R`                             | Runs GSVA using gene expression signatures and a matrix with average gene expressions per cell cluster) |
| `obtains_ORA_for_MatrixColumns.pl`                             | Runs ORA using gene expression signatures and a matrix with average gene expressions per cell cluster) |

**Scripts to run ROC and PR curve analyses**
---
| Script name |  Task(s) |
| -------------------------------------------------------------- |  ------------------------------------------------ |
| `obtains_performance_plots_from_cluster_labelings.pl`             | Compiles results from cell type labeling methods and obtains ROC and PR curves plots and AUC's |
| `obtains_ROC_and_PR_curves_from_matrix_with_gold_standards.R`  | Obtains ROC and PR curve plots, ROC AUC and PR AUC values from a matrix of reference labels in column 2 and predictions in columns 3 to N |

**Scripts to subsample cell type gene expression signatures**
---
| Script name |  Task(s) |
| -------------------------------------------------------------- |  ------------------------------------------------ |
| `obtains_permuted_samples_from_gmt.R`                          | Subsamples genes from gene expression signatures in the form of gene sets |
| `propagates_permuted_gmt_files_to_profile.R`                   | Propagates subsampling from signatures in the form of gene sets to those in the form of gene expression profiles |


**Other scripts**
---
| Script name |  Task(s) |
| -------------------------------------------------------------- |  ------------------------------------------------ |
| `obtains_average_gene_expression_per_cluster.R`                | Obtains a matrix with average gene expressions per cell cluster from scRNA-seq data, and cell cluster assignments) |


How to run the scripts
================
* To see the `help` of R scripts run them like:  <br />
  `Rscript ~/path_to_script/script.R -h`  <br />
  
* To see the help of Perl scripts, make the files executable with  <br />
  `chmod +x script.pl` and run them like:  <br />
  `~/path_to_script/script.pl`  <br />
  
  
Dependencies
================
* Check each script code for dependencies and further documentation.

* To install R script dependencies use `install.package()`,  <br />
  or `install.packages("BiocManager")` and `BiocManager::install()`  <br />
  The full list of R packages is:  <br />
  `"optparse", "vioplot", "GSA", "parallel", "data.table", "GSVA", "qvalue", "precrec", "ROCR", "Seurat", "dplyr", "Rserve", "fisher.test", "p.adjust", "e1071", "preprocessCore", "colorRamps"`  <br />
  Used R version 3.5.1  <br />

* To install Perl script dependencies download `perl_modules` directory from this repository  <br />
  and add it to your `PERL5LIB` environment variable.  <br />
  Other Perl modules required are: `Date::Calc` which can be installed from CPAN  <br />
  Used Perl version 5  <br />

* The following Java scripts are needed: <br />
  CIBERSORT.jar can be obtained from https://cibersort.stanford.edu/download.php  <br />
  gsea-3.0.jar  can be obtained from http://software.broadinstitute.org/gsea/downloads.jsp  <br />
  Used Java version 1.8.0_162 <br />

  
Input Datasets
================
* Three scRNA-seq datasets that can be used as inputs for these scrips, from liver cells (MacParland et al, 2018), peripheral blood mononuclear cells (PBMCs) (Zheng et al, 2017) and retinal neurons (Shekhar et al, 2016), were processed, curated and deposited into Zendo:
https://doi.org/10.5281/zenodo.2575050


Archived code at time of publication
================
http://doi.org/10.5281/zenodo.2583161
    
    
Issues and feature requests
================
Please use this repository templates available at:<br />
[.github/ISSUE_TEMPLATE](https://github.com/jdime/scRNAseq_cell_cluster_labeling/tree/master/.github/ISSUE_TEMPLATE)


Authors
================
Javier Diaz (https://github.com/jdime)
