# universalmotif #

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be
filtered during import based on parameters provided by each motif formats, as
well as during export. Motifs can be export into most major motif formats from
various classes as defined by other Bioconductor packages.

## TODO ##

  - **start using `vapply` instead of `sapply`**

### read_meme ###

  - ~~add compatibility for full MEME format files~~

### various ###

  - complete write_meme
  - read_homer and write_homer
  - meme and homer wrappers (for Linux/Mac/BSD systems only)
  - trim_motifs
  - filter_motifs
  - convert_motifs (including a reverse compliment option via
    `Biostrings::reverseComplement`)
      + `seqLogo` pwm
      + `Biobase` pwm
      + `motifRG` Motif (this one contains a lot of info!)
          * maybe write some methods for this one (e.g. `seqLogo::seqLogo`)
      + `motifStack` pcm and pfm
      + cont.
  - support for other major formats: (as seen on
    [MEME's motif_conversion page](http://meme-suite.org/doc/motif_conversion.html))
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
  - use `Biostrings::PWM` to convert pfm to pwm
  - `motifStack::pcm2pfm`
  - `motifStack::plotMotifLogo(matrixpfm)` --> use this rather than
    `seqLogo::seqLogo`
