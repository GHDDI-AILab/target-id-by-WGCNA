
<!-- README.md is generated from README.Rmd. Please edit that file -->

# targidcn <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

Target identification is an essential first step in drug discovery. This
package implements convenient functions for performing target
identification tasks on gene expression data using the WGCNA method.
(“cn” in the package name stands for “correlation network”.)

## Authors

-   Chen Liang <https://github.com/dzyim>

## Installation

You can install the development version of `targidcn` from
[GitHub](https://github.com/GHDDI-AILab/target-id-by-WGCNA).

``` r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GO.db", "preprocessCore", "impute"))  # Prerequisites for the WGCNA package
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))

if (!require("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("GHDDI-AILab/target-id-by-WGCNA")
```

## Input data formats

The method can be applied to gene expression data generate by RNA-seq or
Mass Spectrometry (MS).

### RNA-seq data

(preprocessing not available so far)

### Labelled MS data

The expression levels of proteins are stored in the columns with the
prefixes `Ratio H/L` and `Ratio H/L normalized`.

### Label-free MS data

The expression levels of proteins are stored in the columns with the
prefixes `LFQ intensity` or `Intensity`.

## Details

R package `targidcn` contains:

-   Functions for loading raw data:
    -   `ReadExperimentalDesign()`: read raw data path, return an object
        of class `ExperimentInfo`
    -   `ReadPhenotypeTable()`: read a phenotype table, return an object
        of class `ExperimentInfo`
    -   `ReadProteinGroups()`: read raw data path, return an object of
        class `ProteinGroups`
-   S3 classes:
    -   `ExperimentList-class`
        -   `ExperimentInfo-class`: for storing sample information
        -   `ExpAssayTable-class`: for storing raw or QCed expression
            data, whose rows correspond to genes
            -   `ProteinGroups-class`: for storing raw or QCed
                proteomics data
        -   `ExpAssayFrame-class`: for storing QCed, scaled, normalized
            expression data, whose rows correspond to samples
            -   `CorrelationNetwork-class`: for storing storing
                expression data with correlation network(s)
-   S3 methods:
    -   `AddPhenotype()`
    -   `Subset()`
    -   `Tidy()`
    -   `QC()`
    -   `Reshape()`: convert an `ExpAssayTable` object to an
        `ExpAssayFrame` object
    -   `LogTransform()`
    -   `Normalize()`
    -   `Histogram()`
    -   `SampleTree()`
    -   `PickThreshold()`
    -   `AddNetwork()`
    -   `ModulePlot()`
    -   `AddConnectivity()`
    -   `GetConnectivity()`
    -   `GetHubGenes()`
    -   `ModuleSignificance()`
    -   `BindModuleSignificance()`
    -   `ModuleTraitHeatmap()`
    -   `GeneSignificance()`
    -   `ModuleMembership()`
    -   `GetRelatedHubGenes()`

## Tutorial

``` r
library(data.table)
library(magrittr)
library(targidcn) %>% suppressMessages()

datadir = system.file("extdata", "MS_label-free", "MaxQuantOutput_50", package = "targidcn")
pheno = ReadPhenotypeTable(file.path(datadir, "phenotype.txt"))
assay = ReadProteinGroups(datadir)
assay
#> An object of class ProteinGroups
#> 
#> 111 Experiment(s): "UCA_1", "UCA_2", "UCA_3", ...
#> 1 Assay(s): "Intensity"
#>  853 features across 111 samples within assay 1.
```

``` r
cn = assay %>% 
  AddPhenotype(pheno) %>% 
  Tidy() %>% 
  QC() %>% 
  Reshape() %>% 
  LogTransform() %>% 
  Normalize() %>% 
  AddNetwork(power = 6) %>% 
  AddConnectivity()
#> Warning in QC.ProteinGroups(.): No column for checking false hits!
#> Warning in QC.ProteinGroups(.): No column for checking unique peptides!
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
#>  ..cutHeight not given, setting it to 0.999  ===>  99% of the (truncated) height range in dendro.
#>  ..done.
#>  mergeCloseModules: Merging modules whose distance is less than 0.15
#>    Calculating new MEs...
```

``` r
cn
#> An object of class CorrelationNetwork
#> 
#> 111 Experiment(s): "UCA_1", "UCA_2", "UCA_3", ...
#> 1 Assay(s): "Intensity"
#>  846 features across 111 samples within assay 1.
#> 
#> Attributes:
#> List of 5
#>  $ phenotype    :'data.frame':   111 obs. of  3 variables:
#>   ..$ Diagnosis: chr [1:111] "UC" "UC" "UC" "UC" ...
#>   ..$ Illness  : int [1:111] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Remission: int [1:111] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ QC           :Classes 'data.table' and 'data.frame':  1 obs. of  6 variables:
#>   ..$ Assay                         : chr "Intensity"
#>   ..$ Raw data                      : int 853
#>   ..$ Remove false hits             : int 853
#>   ..$ With gene names               : int 853
#>   ..$ Unique peptides >= 2          : int 853
#>   ..$ goodGenes, min.fraction >= 0.5: int 846
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ powerEstimate: list()
#>  $ network      :List of 1
#>   ..$ Intensity:List of 10
#>   .. ..$ power           : num 6
#>   .. ..$ MEDissThres     : num 0.15
#>   .. ..$ minModuleSize   : num 30
#>   .. ..$ adjacency       : num [1:846, 1:846] 1.00 1.78e-07 3.47e-04 1.31e-06 2.02e-12 ...
#>   .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. ..$ : chr [1:846] "A1BG" "AGR2" "ANXA2.ANXA2P2" "BSG" ...
#>   .. .. .. ..$ : chr [1:846] "A1BG" "AGR2" "ANXA2.ANXA2P2" "BSG" ...
#>   .. ..$ dissTOM         : num [1:846, 1:846] 0 0.999 0.997 0.999 1 ...
#>   .. ..$ geneTree        :List of 7
#>   .. .. ..$ merge      : int [1:845, 1:2] -17 -352 -23 -18 -318 -597 -431 -493 3 -173 ...
#>   .. .. ..$ height     : num [1:845] 0.62 0.677 0.686 0.687 0.691 ...
#>   .. .. ..$ order      : int [1:846] 796 250 187 598 175 234 628 712 58 227 ...
#>   .. .. ..$ labels     : NULL
#>   .. .. ..$ method     : chr "average"
#>   .. .. ..$ call       : language fastcluster::hclust(d = stats::as.dist(dissTOM), method = "average")
#>   .. .. ..$ dist.method: NULL
#>   .. .. ..- attr(*, "class")= chr "hclust"
#>   .. ..$ moduleEigengenes:'data.frame':  111 obs. of  5 variables:
#>   .. .. ..$ MEturquoise: num [1:111] -0.2971 -0.0318 0.0856 -0.132 0.0428 ...
#>   .. .. ..$ MEblue     : num [1:111] 0.0794 0.063 -0.0601 -0.0131 -0.2603 ...
#>   .. .. ..$ MEbrown    : num [1:111] -0.024 0.076 -0.0939 -0.1834 -0.1365 ...
#>   .. .. ..$ MEyellow   : num [1:111] -0.1034 0.061 -0.0071 -0.0687 -0.1588 ...
#>   .. .. ..$ MEgrey     : num [1:111] -0.1788 -0.1398 -0.0985 -0.1133 -0.1283 ...
#>   .. ..$ moduleColors    : Named chr [1:846] "turquoise" "blue" "yellow" "turquoise" ...
#>   .. .. ..- attr(*, "names")= chr [1:846] "A1BG" "AGR2" "ANXA2.ANXA2P2" "BSG" ...
#>   .. ..$ moduleLabels    : Named num [1:846] 1 2 4 1 0 1 4 1 1 1 ...
#>   .. .. ..- attr(*, "names")= chr [1:846] "A1BG" "AGR2" "ANXA2.ANXA2P2" "BSG" ...
#>   .. ..$ unmergedColors  : Named chr [1:846] "turquoise" "blue" "yellow" "turquoise" ...
#>   .. .. ..- attr(*, "names")= chr [1:846] "A1BG" "AGR2" "ANXA2.ANXA2P2" "BSG" ...
#>  $ connectivity :List of 1
#>   ..$ Intensity:Classes 'data.table' and 'data.frame':   846 obs. of  6 variables:
#>   .. ..$ gene   : chr [1:846] "C3" "HP" "APOA1" "A2M" ...
#>   .. ..$ kTotal : num [1:846] 0.915 1.144 1.013 0.553 0.667 ...
#>   .. ..$ kWithin: num [1:846] 0.424 0.379 0.311 0.176 0.122 ...
#>   .. ..$ kOut   : num [1:846] 0.49 0.765 0.701 0.377 0.545 ...
#>   .. ..$ kDiff  : num [1:846] -0.0659 -0.3865 -0.3899 -0.2008 -0.4232 ...
#>   .. ..$ module : num [1:846] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, ".internal.selfref")=<externalptr>
```

``` r
cn %>% Histogram(preview = TRUE)
#> Warning: Removed 14993 rows containing non-finite values (stat_bin).
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="60%" />

``` r
cn %>% SampleTree(preview = TRUE)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="90%" />

``` r
cn %>% ModulePlot(preview = TRUE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="90%" />

``` r
cn %>% GetHubGenes()
#>      gene         ensembl                                             fullname
#> 1: ATP12A ENSG00000075673 ATPase H+/K+ transporting non-gastric alpha2 subunit
#> 2: ATP1A1 ENSG00000163399           ATPase Na+/K+ transporting subunit alpha 1
#> 3:   ENO1 ENSG00000074800                                            enolase 1
#> 4: COL6A1 ENSG00000142156                       collagen type VI alpha 1 chain
#> 5:    DCN ENSG00000011465                                              decorin
#> 6: COL6A2 ENSG00000142173                       collagen type VI alpha 2 chain
#> 7:    DLD ENSG00000091140                       dihydrolipoamide dehydrogenase
#> 8: SUCLG2 ENSG00000172340        succinate-CoA ligase GDP-forming subunit beta
#>       kTotal   kWithin      kOut      kDiff module
#> 1:  8.493045  5.030402 3.4626427  1.5677597      1
#> 2:  5.720120  4.883286 0.8368346  4.0464514      1
#> 3:  8.331755  6.325000 2.0067549  4.3182450      2
#> 4: 14.027675 11.122084 2.9055901  8.2164943      3
#> 5: 14.676139 10.986816 3.6893238  7.2974918      3
#> 6: 13.141510 10.078089 3.0634212  7.0146680      3
#> 7:  8.740479  4.281549 4.4589306 -0.1773816      4
#> 8:  7.336496  3.912245 3.4242513  0.4879937      4
```

``` r
samples1 = pheno$table[is.na(Remission) | Remission == 0, Experiment]
samples2 = pheno$table[!is.na(Remission), Experiment]
mt1 = ModuleSignificance(cn, samples = samples1, traits = "Illness", prefix = "UC")
mt2 = ModuleSignificance(cn, samples = samples2, traits = "Remission", prefix = "UC")
mt = BindModuleSignificance(mt1, mt2)
```

``` r
mt %>% ModuleTraitHeatmap(preview = TRUE)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="60%" />

``` r
mt %>% GetRelatedHubGenes()
#>      gene         ensembl                                             fullname
#> 1: ATP12A ENSG00000075673 ATPase H+/K+ transporting non-gastric alpha2 subunit
#> 2: ATP1A1 ENSG00000163399           ATPase Na+/K+ transporting subunit alpha 1
#> 3:    DLD ENSG00000091140                       dihydrolipoamide dehydrogenase
#> 4: SUCLG2 ENSG00000172340        succinate-CoA ligase GDP-forming subunit beta
#>      kTotal  kWithin      kOut      kDiff module
#> 1: 8.493045 5.030402 3.4626427  1.5677597      1
#> 2: 5.720120 4.883286 0.8368346  4.0464514      1
#> 3: 8.740479 4.281549 4.4589306 -0.1773816      4
#> 4: 7.336496 3.912245 3.4242513  0.4879937      4
```

## References

**Analysis of oncogenic signaling networks in glioblastoma identifies
ASPM as a molecular target.**  
Horvath S, Zhang B, Carlson M, et al.  
PNAS. 2006;103(46):17402-17407. <doi:10.1073/pnas.0608396103>

**WGCNA: an R package for weighted correlation network analysis.**  
Langfelder P, Horvath S.  
BMC Bioinformatics. 2008;9:559. <doi:10.1186/1471-2105-9-559>

**Structural weakening of the colonic mucus barrier is an early event in
ulcerative colitis pathogenesis.**  
van der Post S, Jabbar KS, Birchenough G, et al.  
Gut. 2019;68(12):2142-2151. <doi:10.1136/gutjnl-2018-317571>

<!-- 
You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>.
--->
