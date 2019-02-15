# corgi: Correlation matrix singular value Gap Inflation

CorGI a gene filter for comparative analysis between two single-cell RNAseq datasets. Conceptually, CorGI selects genes that highlight the _shared heterogenity_ in both datasets while rejecting genes that contributes to _batch effects_.

For example, PCA run on _genes selected by CorGI_ can remove the "batch effect" between mouse and human pre-implantation embryogenesis:


<img src="https://yutongwangumich.github.io/corgi/articles/corgi_files/figure-html/unnamed-chunk-18-1.png" alt="drawing" width = "60%"/>


In contrast, PCA run on _all genes_ results in significant batch effect, even with batch normalization:

<img src="https://yutongwangumich.github.io/corgi/articles/corgi_files/figure-html/unnamed-chunk-17-1.png" alt="drawing" width = "60%"/>


## Installation

In `R`, run

```
library(devtools)
install_github("YutongWangUMich/corgi")
```

## Tutorial

For a quick introduction, check out the vignette: [comparative analysis of pre-implantation embryogenesis between mouse and human](https://yutongwangumich.github.io/corgi/articles/corgi.html).

