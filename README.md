An introduction to Rbwa
================

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Overview](#overview)
-   [Build a reference index with
    `bwa_build_index`](#build-a-reference-index-with-bwa_build_index)
-   [Aligment with `bwa_aln`](#aligment-with-bwa_aln)
    -   [Creating a SAM file](#creating-a-sam-file)
    -   [Creating a SAM file with secondary
        alignments](#creating-a-sam-file-with-secondary-alignments)
-   [Aligment with `bwa_mem`](#aligment-with-bwa_mem)
-   [Session info](#session-info)
-   [References](#references)

Authors: Jean-Philippe Fortin

Date: July 13, 2022

# Introduction

The `Rbwa` package provides an **R** wrapper around the two popular
*BWA* aligners `BWA-backtrack` (Li and Durbin 2009) and `BWA-MEM` (Li
2013).

As mentioned in the BWA manual (see
<http://bio-bwa.sourceforge.net/bwa.shtml>), BWA-backtrack is designed
for short Illumina reads (up to 100bp), while BWA-MEM is more suitable
for longer sequences (70bp to 1Mbp) and supports split alignment.

# Installation

`Rbwa` can be installed from Bioconductor using the following commands
in a fresh R session:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rbwa")
```

# Overview

The two main alignment functions are:

-   BWA-backtrack: `bwa_aln`
-   BWA-MEM: `bwa_mem`

The package also includes the following convenience functions:

-   Genome indexing: `bwa_build_index`
-   Conversion of `bwa_aln` output into SAM output: `bwa_sam`
-   Generation of secondary alignments: `xa2multi`

# Build a reference index with `bwa_build_index`

Both `bwa_aln` and `bwa_mem` require first to create a genome index from
a FASTA file. This is done only once for a given genome. This can be
done using the function `bwa_build_index`.

First, we load `Rbwa`:

``` r
library(Rbwa)
```

In the following example code, we build a BWA index for a small portion
of human chr12 that we stored in a FASTA file located within the `Rbwa`
package. We store the index files in a temporary directory.

``` r
dir <- tempdir()
fasta <- system.file(package="Rbwa",
                     "fasta/chr12.fa")
index_prefix <- file.path(dir, "chr12")
bwa_build_index(fasta,
                index_prefix=index_prefix)
list.files(dir)
```

    ## [1] "chr12.amb" "chr12.ann" "chr12.bwt" "chr12.pac" "chr12.sa"

# Aligment with `bwa_aln`

We now align read sequences stored in the toy example FASTQ file
`fastq/sequences.fastq`, provided in the `Rbwa` package, to our indexed
genome:

``` r
fastq <- system.file(package="Rbwa",
                     "fastq/sequences.fastq")
bwa_aln(index_prefix=index_prefix,
        fastq_files=fastq,
        sai_files=file.path(dir, "output.sai"))
```

Any valid BWA arguments can be passed to the `bwa_aln` function. To see
the complete list of valid arguments, please visit the BWA reference
manual: <http://bio-bwa.sourceforge.net/bwa.shtml>.

For instance, we can specify the maximal edit distance between the query
sequence and the reference genome to be 3 using `n`, as well as the
maximal edit distance in the seed sequence `k` to be 3, where we specify
that the length of the seed sequence is 13 using the argument `l`:

``` r
bwa_aln(index_prefix=index_prefix,
        fastq_files=fastq,
        sai_files=file.path(dir, "output.sai"),
        n=3, k=3, l=13)
```

## Creating a SAM file

The output of `bwa_aln` is an intermediate `sai` file that should be
converted into a `sam` file using the `bwa_sam` function as follows:

``` r
bwa_sam(index_prefix=index_prefix,
        fastq_files=fastq,
        sai_files=file.path(dir, "output.sai"),
        sam_file=file.path(dir, "output.sam"))
```

Let’s read the first few lines of the SAM file:

``` r
strtrim(readLines(file.path(dir, "output.sam")), 65)
```

    ## [1] "@SQ\tSN:chr12\tLN:171330"                                                             
    ## [2] "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1198-dirty\tCL:/Library/Frameworks/R.fram"            
    ## [3] "ACATCAGAAAGAGCGGCAG\t0\tchr12\t170881\t37\t19M\t*\t0\t0\tACATCAGAAAGAGCGGCAG\t~~~~~~~"
    ## [4] "CAACCCAGCCCCCCTCCAA\t0\tchr12\t170801\t37\t19M\t*\t0\t0\tCAACCCAGCCCCCCTCCAA\t~~~~~~~"
    ## [5] "CCTGTGATCCACGGAGGCT\t0\tchr12\t170765\t37\t19M\t*\t0\t0\tCCTGTGATCCACGGAGGCT\t~~~~~~~"
    ## [6] "GCACTGCGGTGAGTGCTGT\t0\tchr12\t170665\t37\t19M\t*\t0\t0\tGCACTGCGGTGAGTGCTGT\t~~~~~~~"
    ## [7] "GCCTTTTACAGTTCGTACT\t0\tchr12\t170820\t37\t19M\t*\t0\t0\tGCCTTTTACAGTTCGTACT\t~~~~~~~"
    ## [8] "GTCATGCCCCCTCAGCCAG\t0\tchr12\t170703\t37\t19M\t*\t0\t0\tGTCATGCCCCCTCAGCCAG\t~~~~~~~"
    ## [9] "TCGGCTCTCACCGTGTCCG\t0\tchr12\t170646\t37\t19M\t*\t0\t0\tTCGGCTCTCACCGTGTCCG\t~~~~~~~"

## Creating a SAM file with secondary alignments

By default, each row of the SAM output corresponds to the best alignment
hit for a given input query sequence. Other alignments (secondary
alignments, or other loci in case of multiple alignments) are stored in
the XA tag.

The function `xa2multi` conveniently extracts the alignments from the XA
tags and represent them as additional rows in the SAM format. This can
be executed as follows:

``` r
xa2multi(file.path(dir, "output.sam"),
         file.path(dir, "output.multi.sam"))
strtrim(readLines(file.path(dir, "output.multi.sam")), 65)
```

    ## [1] "@SQ\tSN:chr12\tLN:171330"                                                             
    ## [2] "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1198-dirty\tCL:/Library/Frameworks/R.fram"            
    ## [3] "ACATCAGAAAGAGCGGCAG\t0\tchr12\t170881\t37\t19M\t*\t0\t0\tACATCAGAAAGAGCGGCAG\t~~~~~~~"
    ## [4] "CAACCCAGCCCCCCTCCAA\t0\tchr12\t170801\t37\t19M\t*\t0\t0\tCAACCCAGCCCCCCTCCAA\t~~~~~~~"
    ## [5] "CCTGTGATCCACGGAGGCT\t0\tchr12\t170765\t37\t19M\t*\t0\t0\tCCTGTGATCCACGGAGGCT\t~~~~~~~"
    ## [6] "GCACTGCGGTGAGTGCTGT\t0\tchr12\t170665\t37\t19M\t*\t0\t0\tGCACTGCGGTGAGTGCTGT\t~~~~~~~"
    ## [7] "GCCTTTTACAGTTCGTACT\t0\tchr12\t170820\t37\t19M\t*\t0\t0\tGCCTTTTACAGTTCGTACT\t~~~~~~~"
    ## [8] "GTCATGCCCCCTCAGCCAG\t0\tchr12\t170703\t37\t19M\t*\t0\t0\tGTCATGCCCCCTCAGCCAG\t~~~~~~~"
    ## [9] "TCGGCTCTCACCGTGTCCG\t0\tchr12\t170646\t37\t19M\t*\t0\t0\tTCGGCTCTCACCGTGTCCG\t~~~~~~~"

# Aligment with `bwa_mem`

The `bwa_mem` function works similar to the `bwa_aln` function, except
that it does not produce intermediate `.sai` files; it outputs a SAM
file directly:

``` r
fastq <- system.file(package="Rbwa",
                     "fastq/sequences.fastq")
bwa_mem(index_prefix=index_prefix,
        fastq_files=fastq,
        sam_file=file.path(dir, "output.sam"))
```

``` r
strtrim(readLines(file.path(dir, "output.sam")), 65)
```

    ## [1] "@SQ\tSN:chr12\tLN:171330"                                                               
    ## [2] "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1198-dirty\tCL:/Library/Frameworks/R.fram"              
    ## [3] "ACATCAGAAAGAGCGGCAG\t4\t*\t0\t0\t*\t*\t0\t0\tACATCAGAAAGAGCGGCAG\t~~~~~~~~~~~~~~~~~~~\t"
    ## [4] "CAACCCAGCCCCCCTCCAA\t4\t*\t0\t0\t*\t*\t0\t0\tCAACCCAGCCCCCCTCCAA\t~~~~~~~~~~~~~~~~~~~\t"
    ## [5] "CCTGTGATCCACGGAGGCT\t4\t*\t0\t0\t*\t*\t0\t0\tCCTGTGATCCACGGAGGCT\t~~~~~~~~~~~~~~~~~~~\t"
    ## [6] "GCACTGCGGTGAGTGCTGT\t4\t*\t0\t0\t*\t*\t0\t0\tGCACTGCGGTGAGTGCTGT\t~~~~~~~~~~~~~~~~~~~\t"
    ## [7] "GCCTTTTACAGTTCGTACT\t4\t*\t0\t0\t*\t*\t0\t0\tGCCTTTTACAGTTCGTACT\t~~~~~~~~~~~~~~~~~~~\t"
    ## [8] "GTCATGCCCCCTCAGCCAG\t4\t*\t0\t0\t*\t*\t0\t0\tGTCATGCCCCCTCAGCCAG\t~~~~~~~~~~~~~~~~~~~\t"
    ## [9] "TCGGCTCTCACCGTGTCCG\t4\t*\t0\t0\t*\t*\t0\t0\tTCGGCTCTCACCGTGTCCG\t~~~~~~~~~~~~~~~~~~~\t"

# Session info

``` r
sessionInfo()
```

    ## R Under development (unstable) (2022-03-21 r81954)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] Rbwa_1.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.2.0  magrittr_2.0.2  fastmap_1.1.0   cli_3.3.0      
    ##  [5] tools_4.2.0     htmltools_0.5.2 rstudioapi_0.13 yaml_2.3.5     
    ##  [9] stringi_1.7.6   rmarkdown_2.13  knitr_1.37      stringr_1.4.0  
    ## [13] xfun_0.30       digest_0.6.29   rlang_1.0.2     evaluate_0.15

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bwa2" class="csl-entry">

Li, Heng. 2013. “Aligning Sequence Reads, Clone Sequences and Assembly
Contigs with BWA-MEM.” *arXiv Preprint arXiv:1303.3997*.

</div>

<div id="ref-bwa1" class="csl-entry">

Li, Heng, and Richard Durbin. 2009. “Fast and Accurate Short Read
Alignment with Burrows–Wheeler Transform.” *Bioinformatics* 25 (14):
1754–60.

</div>

</div>
