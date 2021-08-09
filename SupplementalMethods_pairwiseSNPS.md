"Genomic Stability and Genetic Defense Systems in *Dolosigranulum
pigrum* a Candidate Beneficial Bacterium from the Human Microbiome"
Supplemental Methods
================

# TREE SNP DISTANCE ANALYSIS

## Get the pairwise distance information for the just *D. pigrum* phylogeny

``` r
# read in the nucleotide alignment 
dpi_aln<-read.FASTA("Dpi_concat_core_CDSalignment_taxa28.fa", type="DNA")

# read in the IQTree
dpi_tree<-read.tree("iqtree_CDSnrnogb_modGTRFR3_v2_t28_treefile.tre")
```

Using the R package `harrietr`, we can calculate the evolutionary
pairwise SNP distance in this phylogeny

``` r
dist_df<-dist_long(dpi_aln, tree = dpi_tree)
```

Save the resulting files in a CSV
file

``` r
write.csv(dist_df, "Dpi_SNPdata_unrootedtree.csv")
```

## Get the pairwise distance information for the tree with *A. otitis* ATCC 51267 as an outgroup

``` r
# read in the nucleotide alignment 
dpi_rtaln<-read.FASTA("Dpi_concat_core_CDSalignment_outgroupAO_taxa29.fa", type="DNA")
# read in the IQTree
dpi_rttree<-read.tree("iqtree_CDSrtAOghnogb_t29_contree.tre")

#calculate the pairwise SNP distance using harrietr
dist_df3<-dist_long(dpi_rtaln, tree = dpi_rttree)
```

Save the resulting files in a CSV
file

``` r
write.csv(dist_df3, "Dpi_SNPdata_rootedtree.csv")
```

## Combining the results from the roote and unrooted tree into one file that will become our **Table S1a**

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
```

Save the resulting files in a CSV file

``` r
write.csv(comb_distf, "Dpi_SNPdata_bothstrees.csv")
```
