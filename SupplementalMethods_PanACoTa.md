"Genomic Stability and Genetic Defense Systems in *Dolosigranulum
pigrum* a Candidate Beneficial Bacterium from the Human Microbiome"
Supplemental Methods
================

# PanACoTa ANALYSIS

[PanACoTa, v1.2.0](https://github.com/gem-pasteur/PanACoTA) ([Perrin and
Rocha 2021](#ref-10.1093/nargab/lqaa106)) was installed in a Python
environment called `PanACoTa`.

## Mash Analysis

We used `PanACoTa` in order to calculate `Mash` distances ([Ondov et al.
2016](#ref-Ondov2016)) between each pair of genomes:

``` bash
#conda activate PanACoTa
PanACoTA prepare --norefseq -o analysis_PanACoTa -d GENOMES/renamed
```

# <u>REFERENCES</u>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Ondov2016" class="csl-entry">

Ondov, Brian D., Todd J. Treangen, Páll Melsted, Adam B. Mallonee,
Nicholas H. Bergman, Sergey Koren, and Adam M. Phillippy. 2016. “Mash:
Fast Genome and Metagenome Distance Estimation Using MinHash.” *Genome
Biology* 17 (1). <https://doi.org/10.1186/s13059-016-0997-x>.

</div>

<div id="ref-10.1093/nargab/lqaa106" class="csl-entry">

Perrin, Amandine, and Eduardo P C Rocha. 2021. “<span
class="nocase">PanACoTA: a modular tool for massive microbial
comparative genomics</span>.” *NAR Genomics and Bioinformatics* 3 (1).
<https://doi.org/10.1093/nargab/lqaa106>.

</div>

</div>
