"Genomic Stability and Genetic Defense Systems in *Dolosigranulum
pigrum* a Candidate Beneficial Bacterium from the Human Microbiome"
Supplemental Methods
================

# SELECTED MGEs FROM THE ANVI’O & PPanGGOLiN ANALYSIS

## Anvi’o 7 Annotation search

We import the output of `anvi-summarize` and select the most relevant
variables for this analysis:

``` r
DpigPangenome <-  read_delim("analysis_Anvio7/Pangenomic_Results_Dpig/Dpig-PAN-SUMMARY/PAN_DPIG_prokka_gene_clusters_summary.txt.gz", "\t", escape_double = FALSE, trim_ws = TRUE)
DpigPangenome <- DpigPangenome %>%
  select(-unique_id, -aa_sequence, -SCG, -functional_homogeneity_index, -geometric_homogeneity_index, -combined_homogeneity_index)
```

We search for individual genes annotated as **Transposases**/**Retrons**
in either the Prokka, COG20, Pfam or KOfam annotations:

``` r
Retron_Genes <- DpigPangenome %>%
  filter(str_detect(COG20_FUNCTION, "Retron") | str_detect(`Prokka:Prodigal`, "retron") | str_detect(Pfam, "Retron") | str_detect(KOfam, "retron"))
Intron_Genes <- DpigPangenome %>%
  filter(str_detect(COG20_FUNCTION, "Intron") | str_detect(`Prokka:Prodigal`, "intron") | str_detect(Pfam, "Intron") | str_detect(KOfam, "intron"))
Transposase_Genes <- DpigPangenome %>%
  filter(str_detect(COG20_FUNCTION, "Transposase") | str_detect(`Prokka:Prodigal`, "transposase") | str_detect(Pfam, "Transposase")| str_detect(KOfam, "transposase"))
```

There are 153 individual genes annotated as **Retrons**, 3 as
**Introns** and 299 as **Transposases**.

We identified the GCs that contain those individual gene annotations:

``` r
Retron_GCs <- Retron_Genes %>%
  group_by(gene_cluster_id) %>%
  summarize(n = n())
Intron_GCs <- Intron_Genes %>%
  group_by(gene_cluster_id) %>%
  summarize(n = n())
Transposases_GCs <- Transposase_Genes %>%
  group_by(gene_cluster_id) %>%
  summarize(n = n())
```

There are 3 GCs annotated as **Retrons**, 1 as **Introns** and 23 as
**Transposases**.

### Generating Anvi’o fasta files

From Anvi’o we export the protein Sequences for the GCs with individual
genes identified as **Introns**/**Retrons** and as **Transposases**:

For the **Introns**/**Retrons** we focused only on the main GC:
GC\_00000001:

``` bash
#conda activate anvio-7
anvi-get-sequences-for-gene-clusters -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db \
                                     -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db \
                                     -o "analysis_MGEs/SequencesAnvio7/GC_00000001_Retron.faa" \
                                     --gene-cluster-id GC_00000001
```

For the **Transposases** we first write a file with the GC numbers for
the 23 identified GCs:

``` r
write.csv(Transposases_GCs$gene_cluster_id, "analysis_MGEs/SequencesAnvio7/IDs_Transposases.txt", row.names = FALSE)
```

(File cleaned up on text editor to removed header and "")

``` bash
#conda activate anvio-7
anvi-get-sequences-for-gene-clusters -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db \
                                     -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db \
                                     -o "analysis_MGEs/SequencesAnvio7/GCsTransposases.faa" \
                                     --gene-cluster-ids-file "analysis_MGEs/SequencesAnvio7/IDs_Transposases.txt"
```

### Selecting representative sequences

For each GC alignments were visually inspected in AliView and
full-length representative sequences selected for PFam search. Using the
[PFam batch sequence
search](http://pfam.xfam.org/search#tabview=tab1)/[HMMER
website](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) ([Mistry et
al. 2020](#ref-mistry2020); [Potter et al. 2018](#ref-Potter2018)) we
classified the initial 23 putative **Tranposase** GSs as:

-   **Real Transposases:** GC\_00000003, GC\_00000040, GC\_00000055,
    GC\_00001693, GC\_00002092, GC\_00002210, GC\_00002310 and
    GC\_00002501.
-   **Integrases (rve domain):** GC\_00000028, GC\_00000085,
    GC\_00001701, GC\_00001775 and GC\_00002348.
-   **Other/Partial:** GC\_00000008, GC\_00001669, GC\_00001787,
    GC\_00002105, GC\_00002382, GC\_00002430, GC\_00002460,
    GC\_00002491, GC\_00002679 and GC\_00002805.

## BOFFO/Clinker: UPDATE & ADD REFERENCES HERE

The following files were created with the selected sequences:

-   SelectedAnvio7\_Intron.faa: Representative sequence for the
    GC\_00000001 cluster.
-   SelectedAnvio7\_Real\_Transposases.faa: GCs initially identified
    with the word Transposase on the annotation search and with complete
    (80% coverage or more) PFam **Transposase** domains.
-   SelectedAnvio7\_Integrases\_rve.faa: GCs initially identified with
    the word Transposase on the annotation search but with complete (80%
    coverage or more) PFam **rve** domains.

These files were analyzed with BOFFO with minimum percent identity 85%
and minimum coverage 80%.

Output is in AllMGE\_BOFFO\_80cov85id: You will see a folder for each of
the multi-fasta’s sent to him (“Selected\_Real\_Transposases\_v2” for
example). Within those folders there are subfolders for each protein
(“GC\_00000084,” etc.) with all of the BOFFO outputs. The clinker
results are located in a folder called “html” and if you do not see this
folder it means that only one genome hit was found. Since clinker
requires a minimum of two, there is no output.

Outputs concatenated by group:

``` bash
cd ..
cat AllMGE_BOFFO_80cov85id/Selected_Intron.faa-output/*/*.csv.gz  > SelectedGCs_Anvio/BOFFO/Intron.txt.gz
cat AllMGE_BOFFO_80cov85id/Selected_Real_Transposases_v2-output/*/*.csv.gz  > SelectedGCs_Anvio/BOFFO/Real_Transposases.txt.gz
```

Files unzipped and format cleaned in excel and imported into R:

``` r
Intron <- read_excel("BOFFO/Intron.xlsx")
Real_Transposases <- read_excel("BOFFO/Real_Transposases.xlsx")
ALL <- read_excel("BOFFO/ALL.xlsx")
```

### Final Tables

Summary matrix table with a Genome in each row and each GC listed in
each column.

``` r
MatrixBOFFO_ALL <- ALL %>%
  group_by(genome_name, gene_name) %>%
  summarize(n = n())
MatrixBOFFO_ALL <- spread(MatrixBOFFO_ALL, gene_name, n)
MatrixBOFFO_ALL <- MatrixBOFFO_ALL %>% remove_rownames %>% column_to_rownames(var="genome_name")
```

Totals/Stats by column/row:

``` r
GCTotal <- colSums(MatrixBOFFO_ALL, na.rm=TRUE)
GCMean <- round(colMeans(MatrixBOFFO_ALL, na.rm=TRUE), 2)
GCMedian <- round(colMedians(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2)
GCVar <- round(colVars(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) #Sample variance
GCSD <- round(colSds(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) #Sample standard deviation
GCMAD <- round(colMads(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) #Median absolute deviation
GCMin <- round(colMins(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) 
GCMax <- round(colMaxs(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) 

MatrixBOFFO_ALL["GCTotal" ,] <- GCTotal
MatrixBOFFO_ALL["GCMean" ,] <- GCMean
MatrixBOFFO_ALL["GCMedian" ,] <- GCMedian
MatrixBOFFO_ALL["GCVar" ,] <- GCVar
MatrixBOFFO_ALL["GCSD" ,] <- GCSD
MatrixBOFFO_ALL["GCMAD" ,] <- GCMAD
MatrixBOFFO_ALL["GCMin" ,] <- GCMin
MatrixBOFFO_ALL["GCMax" ,] <- GCMax

MatrixBOFFO_ALL$Total <- rowSums(MatrixBOFFO_ALL, na.rm=TRUE)
```

``` r
write.csv(MatrixBOFFO_ALL, "Stats_MatrixBOFFO.csv", row.names = TRUE)
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-mistry2020" class="csl-entry">

Mistry, Jaina, Sara Chuguransky, Lowri Williams, Matloob Qureshi,
Gustavo A Salazar, Erik L L Sonnhammer, Silvio C E Tosatto, et al. 2020.
“Pfam: The Protein Families Database in 2021.” *Nucleic Acids Research*
49 (D1): D412–19. <https://doi.org/10.1093/nar/gkaa913>.

</div>

<div id="ref-Potter2018" class="csl-entry">

Potter, Simon C, Aurélien Luciani, Sean R Eddy, Youngmi Park, Rodrigo
Lopez, and Robert D Finn. 2018. “HMMER Web Server: 2018 Update.”
*Nucleic Acids Research* 46 (W1): W200–204.
<https://doi.org/10.1093/nar/gky448>.

</div>

</div>
