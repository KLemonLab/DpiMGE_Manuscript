"Genomic Stability and Genetic Defense Systems in *Dolosigranulum
pigrum* a Candidate Beneficial Bacterium from the Human Microbiome"
Supplemental Methods
================

# TREE SNP DISTANCE ANALYSIS

## Pairwise distance information for just the *D. pigrum* phylogeny

``` r
# read in the nucleotide alignment 
dpi_aln<-read.FASTA("analysis_PhylogeneticDistances/SNPs/Dpi_concat_core_CDSalignment_taxa28.fa", type="DNA")

# read in the IQTree
dpi_tree<-read.tree("analysis_PhylogeneticDistances/SNPs/iqtree_CDSnrnogb_modGTRFR3_v2_t28_contree.tre")
```

Using the R package
[`harrietr`](https://cran.r-project.org/web/packages/harrietr/README.html),
we can calculate the evolutionary pairwise SNP distances in this
phylogeny

``` r
dist_df<-dist_long(dpi_aln, tree = dpi_tree)
kable(head(dist_df))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
iso1
</th>
<th style="text-align:left;">
iso2
</th>
<th style="text-align:right;">
dist
</th>
<th style="text-align:right;">
evol\_dist
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1914
</td>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:right;">
21506
</td>
<td style="text-align:right;">
0.0344338
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1931\_CDC4294-98
</td>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:right;">
21576
</td>
<td style="text-align:right;">
0.0346119
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL3077
</td>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:right;">
21229
</td>
<td style="text-align:right;">
0.0345185
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1932\_CDC4420-98
</td>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:right;">
22348
</td>
<td style="text-align:right;">
0.0360584
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1938\_CDC4791-99
</td>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:right;">
23149
</td>
<td style="text-align:right;">
0.0374657
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL3050
</td>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:right;">
23099
</td>
<td style="text-align:right;">
0.0388156
</td>
</tr>
</tbody>
</table>

Save the resulting files in a CSV file

``` r
write.csv(dist_df, "analysis_PhylogeneticDistances/SNPs/Dpi_SNPdata_unrootedtree.csv")
```

## Pairwise distance information for the tree with *A. otitis* ATCC 51267 as an outgroup

``` r
# read in the nucleotide alignment 
dpi_rtaln<-read.FASTA("analysis_PhylogeneticDistances/SNPs/Dpi_concat_core_CDSalignment_outgroupAO_taxa29.fa", type="DNA")
# read in the IQTree
dpi_rttree<-read.tree("analysis_PhylogeneticDistances/SNPs/iqtree_CDSrtAOghnogb_t29_contree.tre")

#calculate the pairwise SNP distance using harrietr
dist_df3<-dist_long(dpi_rtaln, tree = dpi_rttree)
kable(head(dist_df3))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
iso1
</th>
<th style="text-align:left;">
iso2
</th>
<th style="text-align:right;">
dist
</th>
<th style="text-align:right;">
evol\_dist
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276773
</td>
<td style="text-align:right;">
1.780329
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1930\_CDC2949-98
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276726
</td>
<td style="text-align:right;">
1.779906
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1934\_CDC4709-98
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276719
</td>
<td style="text-align:right;">
1.779341
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1939\_CDC4792-99
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276654
</td>
<td style="text-align:right;">
1.780876
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL3250
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276762
</td>
<td style="text-align:right;">
1.780435
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1932\_CDC4420-98
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276682
</td>
<td style="text-align:right;">
1.780915
</td>
</tr>
</tbody>
</table>

Save the resulting files in a CSV file

``` r
write.csv(dist_df3, "analysis_PhylogeneticDistances/SNPs/Dpi_SNPdata_rootedtree.csv")
```

## Combining the results from the rooted and unrooted tree into one file that will become our **Table S1a**

``` r
comb_dist<-full_join(dist_df3,dist_df, by = c("iso1","iso2"))

orphans<-comb_dist %>%
  filter(is.na(dist.x)) %>%
  select(iso2, iso1, dist.y,evol_dist.y)

colnames(orphans)[colnames(orphans)=="iso2"] <- "iso1.x"
colnames(orphans)[colnames(orphans)=="iso1"] <- "iso2"
colnames(orphans)[colnames(orphans)=="iso1.x"] <- "iso1"

comb_dist2<-left_join(dist_df3,dist_df,by = c("iso1","iso2"))
comb_distf<-left_join(comb_dist2,orphans,by = c("iso1","iso2"))

kable(head(comb_distf))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
iso1
</th>
<th style="text-align:left;">
iso2
</th>
<th style="text-align:right;">
dist.x
</th>
<th style="text-align:right;">
evol\_dist.x
</th>
<th style="text-align:right;">
dist.y.x
</th>
<th style="text-align:right;">
evol\_dist.y.x
</th>
<th style="text-align:right;">
dist.y.y
</th>
<th style="text-align:right;">
evol\_dist.y.y
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
D.pigrum\_ATCC\_51524
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276773
</td>
<td style="text-align:right;">
1.780329
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1930\_CDC2949-98
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276726
</td>
<td style="text-align:right;">
1.779906
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1934\_CDC4709-98
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276719
</td>
<td style="text-align:right;">
1.779341
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1939\_CDC4792-99
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276654
</td>
<td style="text-align:right;">
1.780876
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL3250
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276762
</td>
<td style="text-align:right;">
1.780435
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
D.pigrum\_KPL1932\_CDC4420-98
</td>
<td style="text-align:left;">
A.otitis\_ATCC\_51267
</td>
<td style="text-align:right;">
276682
</td>
<td style="text-align:right;">
1.780915
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

Save the resulting files in a CSV file

``` r
write.csv(comb_distf, "analysis_PhylogeneticDistances/SNPs/Dpi_SNPdata_bothstrees.csv")
```

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
