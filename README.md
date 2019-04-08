# corgi

CORGI (Curve-fit Objective Ranking of Gene Importance) is a gene filter for integrating two single-cell RNAseq trajectory datasets. Conceptually, CORGI selects genes that highlight the _shared heterogenity_ in both datasets while rejecting genes that contributes to _batch effects_.

For example, MDS run on genes selected by CORGI can remove the "batch effect" between mouse and human pre-implantation embryogenesis:


<img src="https://yutongwangumich.github.io/corgi/articles/corgi_files/figure-html/comparison-1.png" alt="drawing" width = "60%"/>

## Installation

In `R`, run

```
library(devtools)
install_github("YutongWangUMich/corgi")
```

## Tutorial

For a quick introduction, check out the vignette: [integrating pre-implantation embryogenesis between mouse and human](https://yutongwangumich.github.io/corgi/articles/corgi.html).

