---
title: 'universalmotif: An R package for biological motif analysis'
tags:
- R
- transcription factor binding sites
- biological sequence motifs
- DNA
authors:
- name: Benjamin Jean-Marie Tremblay
  orcid: 0000-0002-7441-2951
  Corresponding: true
  affiliation: 1
affiliations:
- name: Independent Researcher, Spain
  index: 1
date: "18 August 2024"
bibliography: paper.bib
---

# Summary

Sequence motifs are an important concept in molecular biology, as specific repeating patterns in DNA, RNA and proteins form the basis of biological regulation. Identifying and characterizing these motifs is therefore a important part of studying various aspects of cellular processes, such as gene regulation, transcript stability, and protein function. Many programs have been developed over the years to tackle these tasks, though their interoperability remains poor. The universalmotif package has two main goals: to serve as a go-between for most common biological motif programs and Bioconductor packages used by the research community, and to provide a robust set of tools for basic motif analysis and manipulation in R. Tools for motif and sequence manipulation, scanning, enrichment, comparison, shuffling and P-value computation are included.

# Installation

The universalmotif project including its extensive documentation are hosted on [Bioconductor](https://bioconductor.org/packages/universalmotif/), with pre-built binaries available for macOS and Windows (and installation from source available for all platforms). Installation takes place from within R using the BiocManager package, which itself can be installed from [CRAN](https://CRAN.R-project.org/package=BiocManager):

```
install.packages("BiocManager")
BiocManager::install("universalmotif")
```

# Statement of need

Identifying and characterizing biological sequence motifs is an important task in the field of molecular biology, especially for the understanding of gene, transcript and protein regulation. Over the year many programs have been created to tackle this, such as the hugely popular MEME suite [@bailey94] and the TFBSTools R/Bioconductor package [@tan16], as well as various curated databases such as JASPAR [@rauluseviciute24] containing thousands of published motifs maintained by the community. While collectively these efforts have been responsible for significant advances in the field, the proliferation of different formats and the resulting poor interoperability hinders their use. To solve this problem, the universalmotif R/Bioconductor package allows for the import and export of a large number of commonly used motif formats, including CIS-BP [@weirauch14], HOMER [@heinz10], JASPAR [@rauluseviciute24], HOCOMOCO [@vorontsov24], TRANSFAC [@wingender96], UniPROBE [@hume15], and any additional simple formats via the `read_matrix()` and `write_matrix()` functions. Furthermore, existing R/Bioconductor packages providing their own motif classes can be converted to and from the universal `universalmotif` motif class via the `convert_motifs()` function, which include the popular TFBSTools [@tan16] and motifStack [@ou18] packages, among others. While other R/Bioconductor motif-related packages often provide functions to import external motif formats (such as `importMatrix()` from motifStack), their other functions still typically cannot be used with motif classes from different packages. By allowing for all R/Bioconductor packages to interoperate via the `universalmotif` class, the universalmotif package provides a way for users to pick and choose functions from all available R/Bioconductor packages. Various other projects have now made use of this extensive compatibility and flexible motif class, such as memes [@nystrom21], CollecTRI [@mullerdott23], circRNAprofiler [@aufiero20], ASTK [@huang24], and the standalone RSAT matrix-clustering tool [@castromondragon17]. The universalmotif project has been continuously developed over six years and will continue to add more formats and classes when requested by the community.

The universalmotif package also provides a suite of functions for working with motifs and biological sequences, giving researchers the ability to perform most motif-related tasks from within R (which has been embraced by a large number of molecular biologists as the programming environment of choice for bioinformatic analyses). These include functions for manipulating motifs themselves (`create_motif()`, `convert_type()`, `filter_motifs()`, `motif_rc()`, `switch_alph()`, `trim_motifs()`), comparison and merging (`compare_motifs()`, `merge_motifs()`, `merge_similar()`), plotting (`view_motifs()`), motif P-values (`motif_pvalue()`), and sequence scanning, enrichment, and manipulation (`scan_sequences()`, `enrich_motifs()`, `create_sequences()`, `get_bkg()`, `shuffle_sequences()`, `sequence_complexity()`). Many additional utilities are included, all of which are extensively documented. No other R/Bioconductor motif-related package offers such a large set of functions for working with motifs, though for packages which contain specialized methods not provided by the universalmotif package (such as the advanced motif plotting of the motifStack package), they can be used seamlessly alongside the universalmotif package. These functions from the universalmotif package have now seen widespread use by researchers over the past several years, as evidenced by their appearance in a wide range of journals (to name some prominent examples: @jores21; @mikl22; @zeng22; @li23; @meeuse23; @najle23; @gao24; @gibson24; @hawkins24; @hoge24). Future developments of the universalmotif project are aimed to increase the available functionality of the package.

# Acknowledgements

This software has been developed and maintained without any financial support. I am grateful to present and past mentors for encouraging me to pursue this project, including Barbara Moffatt, Andrew Doxey, and Julia Qüesta.

# References


