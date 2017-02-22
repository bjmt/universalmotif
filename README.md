# universalmotif #

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be
filtered during import based on parameters provided by each motif formats, as
well as during export. Motifs can be export into most major motif formats from
various classes as defined by other Bioconductor packages.

## TODO ##

### read_meme ###

  - add compatibility for full MEME format files

### various ###

  - complete write_meme
  - read_homer and write_homer
  - meme and homer wrappers
  - support for other major formats:
      + beeml
      + chen
      + elm
      + iupac
      + jaspar (old, 2014, 2016, sites, CM)
      + matrix
      + nmica
      + priority
      + rna
      + scpd
      + sites
      + taipale
      + tamo
      + transfac (including transfac-like)
      + uniprobe
