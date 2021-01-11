*Dolosigranulum pigrum* COG ANALYSIS
================

In order to use the desired Python environment we need to define the
path in the .Rprofile for this project.

``` r
file.edit(".Rprofile")
```

In the opened .Rprofile file paste this: Sys.setenv(RETICULATE\_PYTHON =
“/Users/isabelfe/opt/anaconda3/envs/anvio-6.2/bin/python”)

## Data Import

We import the output of `anvi-summarize` and select the most relevant
variables for the functional analysis:

``` r
DpigPangenome <-  read_delim("PAN_DPIG_prokka_gene_clusters_summary.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
DpigPangenome <- DpigPangenome %>%
  select(unique_id, gene_cluster_id, bin_name, genome_name, num_genomes_gene_cluster_has_hits, num_genes_in_gene_cluster, `Prokka:Prodigal_ACC`, `Prokka:Prodigal`, COG_CATEGORY, COG_FUNCTION, COG_FUNCTION_ACC)
```

In the new variable “accessory\_vs\_core” we define “Soft/Core” as
“MC\_core”+“SC\_core”+“soft\_core” and “accessory” as “shell”+“cloud”:

``` r
DpigPangenome <- DpigPangenome %>%
  mutate(accessory_vs_core=ifelse(bin_name=="MC_Core"|bin_name=="SC_Core"|bin_name=="SoftCore","Core","Accessory"))
```

## Transposase/Retron search

We search for individual genes annotated as Transposases/Retrons in
either the COG or Prokka annotations:

``` r
Transposase_Genes <- DpigPangenome %>%
  filter(str_detect(COG_FUNCTION, "Transposase") | str_detect(`Prokka:Prodigal`, "transposase"))
Retron_Genes <- DpigPangenome %>%
  filter(str_detect(COG_FUNCTION, "Retron") | str_detect(`Prokka:Prodigal`, "retron"))
```

There are 278 individual genes annotated in either the COG or Prokka
annotations as Transposases and 153 as Retrons. However, they belong to
GCs that contain some extra individual genes that lack annotation, but
should be considered also, since they belong to the same cluster.

The following datasets include all the members of the GCs identified in
the previous step:

``` r
Transposase_Genes_FullGCs <- DpigPangenome[DpigPangenome$gene_cluster_id %in% Transposase_Genes$gene_cluster_id,]
Retron_Genes_FullGCs <- DpigPangenome[DpigPangenome$gene_cluster_id %in% Retron_Genes$gene_cluster_id,]
```

Therefore, taking into account all the members of the selected GCs there
are 356 individual genes identified as Transposases and 153 as Retrons.

### Tables

#### Transposase Tables

Summary matrix table with a Genome in each row and each transposase GC
listed in each column.

``` r
Transposases_Matrix <- Transposase_Genes_FullGCs %>%
  group_by(genome_name, gene_cluster_id) %>%
  summarize(n = n())
Transposases_Matrix <- spread(Transposases_Matrix, gene_cluster_id, n)
Transposases_Matrix <- Transposases_Matrix %>% remove_rownames %>% column_to_rownames(var="genome_name")
```

Totals/Stats by column/row:

``` r
GCTotal <- colSums(Transposases_Matrix, na.rm=TRUE)
GCMean <- round(colMeans(Transposases_Matrix, na.rm=TRUE), 2)
GCMedian <- round(colMedians(as.matrix(Transposases_Matrix), na.rm=TRUE), 2)
GCVar <- round(colVars(as.matrix(Transposases_Matrix), na.rm=TRUE), 2) #Sample variance
GCSD <- round(colSds(as.matrix(Transposases_Matrix), na.rm=TRUE), 2) #Sample standard deviation
GCMAD <- round(colMads(as.matrix(Transposases_Matrix), na.rm=TRUE), 2) #Median absolute deviation
GCMin <- round(colMins(as.matrix(Transposases_Matrix), na.rm=TRUE), 2) 
GCMax <- round(colMaxs(as.matrix(Transposases_Matrix), na.rm=TRUE), 2) 

Transposases_Matrix["GCTotal" ,] <- GCTotal
Transposases_Matrix["GCMean" ,] <- GCMean
Transposases_Matrix["GCMedian" ,] <- GCMedian
Transposases_Matrix["GCVar" ,] <- GCVar
Transposases_Matrix["GCSD" ,] <- GCSD
Transposases_Matrix["GCMAD" ,] <- GCMAD
Transposases_Matrix["GCMin" ,] <- GCMin
Transposases_Matrix["GCMax" ,] <- GCMax

Transposases_Matrix$Total <- rowSums(Transposases_Matrix, na.rm=TRUE)
```

``` r
Transposases_Names <- Transposase_Genes_FullGCs %>%
  group_by(gene_cluster_id, `Prokka:Prodigal`) %>%
  summarize(n = n())
```

    ## `summarise()` regrouping output by 'gene_cluster_id' (override with `.groups` argument)

``` r
write.csv(Transposases_Matrix, "Stats_Transposases.csv", row.names = TRUE)
write.csv(Transposases_Names, "Transposases_Names.csv", row.names = TRUE)
```

#### Retron Tables

We focus on the main retron, GC number GC\_00000001:

``` r
Retron_GC_00000001 <- DpigPangenome[DpigPangenome$gene_cluster_id =="GC_00000001",]

Retron_GC_00000001_byGenome <- Retron_GC_00000001 %>%
  group_by(genome_name) %>%
  summarize(GC_00000001 = n())

Retron_GC_00000001_byGenome <- Retron_GC_00000001_byGenome %>% remove_rownames %>% column_to_rownames(var="genome_name")
```

Totals/Stats:

``` r
GCTotal <- colSums(Retron_GC_00000001_byGenome, na.rm=TRUE)
GCMean <- round(colMeans(Retron_GC_00000001_byGenome, na.rm=TRUE), 2)
GCMedian <- round(colMedians(as.matrix(Retron_GC_00000001_byGenome), na.rm=TRUE), 2)
GCVar <- round(colVars(as.matrix(Retron_GC_00000001_byGenome), na.rm=TRUE), 2) #Sample variance
GCSD <- round(colSds(as.matrix(Retron_GC_00000001_byGenome), na.rm=TRUE), 2) #Sample standard deviation
GCMAD <- round(colMads(as.matrix(Retron_GC_00000001_byGenome), na.rm=TRUE), 2) #Median absolute deviation
GCMin <- round(colMins(as.matrix(Retron_GC_00000001_byGenome), na.rm=TRUE), 2) 
GCMax <- round(colMaxs(as.matrix(Retron_GC_00000001_byGenome), na.rm=TRUE), 2) 

Retron_GC_00000001_byGenome["GCTotal" ,] <- GCTotal
Retron_GC_00000001_byGenome["GCMean" ,] <- GCMean
Retron_GC_00000001_byGenome["GCMedian" ,] <- GCMedian
Retron_GC_00000001_byGenome["GCVar" ,] <- GCVar
Retron_GC_00000001_byGenome["GCSD" ,] <- GCSD
Retron_GC_00000001_byGenome["GCMAD" ,] <- GCMAD
Retron_GC_00000001_byGenome["GCMin" ,] <- GCMin
Retron_GC_00000001_byGenome["GCMax" ,] <- GCMax
```

``` r
write.csv(Retron_GC_00000001_byGenome, "Stats_Retron.csv", row.names = TRUE)
```

### Fasta Files

DNA and Protein Sequences for the Anvio GCs identified as
Transposases/Retrons:

``` bash
cd ..
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     -o "TransposasesRetrons/GCsTransposases.faa" \
                                     --gene-cluster-ids-file "TransposasesRetrons/TransposasesIDs.txt"
                                     
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     --report-DNA-sequences \
                                     -o "TransposasesRetrons/GCsTransposases.fsa" \
                                     --gene-cluster-ids-file "TransposasesRetrons/TransposasesIDs.txt"
```

``` bash
cd ..
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     -o "TransposasesRetrons/GC_00000003_TransposaseISL3.faa" \
                                     --gene-cluster-id GC_00000003
                                     
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     --report-DNA-sequences \
                                     -o "TransposasesRetrons/GC_00000003_TransposaseISL3.fsa" \
                                     --gene-cluster-id GC_00000003
```

``` bash
cd ..
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     -o "TransposasesRetrons/GC_00000008_TransposaseIS3.faa" \
                                     --gene-cluster-id GC_00000008
                                     
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     --report-DNA-sequences \
                                     -o "TransposasesRetrons/GC_00000008_TransposaseIS3.fsa" \
                                     --gene-cluster-id GC_00000008
```

``` bash
cd ..
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     -o "TransposasesRetrons/GC_00000001_Retron.faa" \
                                     --gene-cluster-id GC_00000001
                                     
anvi-get-sequences-for-gene-clusters -g Dpig-GENOMES.db \
                                     -p PAN_DPIG_prokka-PAN.db \
                                     --report-DNA-sequences \
                                     -o "TransposasesRetrons/GC_00000001_Retron.fsa" \
                                     --gene-cluster-id GC_00000001
```
