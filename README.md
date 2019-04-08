# corgi: a gene filter for comparative analysis

CORGI a gene filter for integrating two single-cell RNAseq trajectory datasets. Conceptually, CorGI selects genes that highlight the _shared heterogenity_ in both datasets while rejecting genes that contributes to _batch effects_.

For example, PCA run on genes selected by CorGI can remove the "batch effect" between mouse and human pre-implantation embryogenesis:


<img src="https://yutongwangumich.github.io/corgi/articles/corgi_files/figure-html/unnamed-chunk-16-1.png" alt="drawing" width = "60%"/>

## Installation

In `R`, run

```
library(devtools)
install_github("YutongWangUMich/corgi")
```

## Tutorial

For a quick introduction, check out the vignette: [comparative analysis of pre-implantation embryogenesis between mouse and human](https://yutongwangumich.github.io/corgi/articles/corgi.html).

