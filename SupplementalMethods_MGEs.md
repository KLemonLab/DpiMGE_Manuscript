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
in either the Prokka ([Seemann 2014](#ref-Seemann2014)), COG20 ([Tatusov
1997](#ref-Tatusov1997); [Galperin et al. 2020](#ref-Galperin2021)),
Pfam ([Mistry et al. 2020](#ref-mistry2020)) or KOfam ([Minoru Kanehisa
et al. 2015](#ref-kanehisa2015); [M. Kanehisa 2000](#ref-kanehisa2000))
annotations:

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

## BOFFO/Clinker

The following files were created with the selected sequences and
analyzed with the Bacterial Operon Finder for Functional Organization,
aka [BOFFO](https://github.com/FredHutch/boffo) to identify the gene
neighborhoods in which the selected genes were located across all 28 D.
pigrum genomes

-   **SelectedAnvio7\_Intron.faa**: Representative sequence for the
    GC\_00000001 cluster.
-   **SelectedAnvio7\_Real\_Transposases.faa**: GCs initially identified
    with the word Transposase on the annotation search and with complete
    (80% coverage or more) PFam **Transposase** domains
-   **SelectedAnvio7\_Integrases\_rve.faa**: GCs initially identified
    with the word Transposase on the annotation search but with complete
    (80% coverage or more) PFam **rve** domains:

The groups of genes identified with BOFFO at minimum percent identity
85% and minimum coverage 80% were visualized using
[clinker](https://github.com/gamcil/clinker) ([Gilchrist et al.
2021](#ref-Gilchrist2021)). The following links display the whole
clinker output:

-   **Intron:**
    [GC\_00000001](https://htmlpreview.github.io/?https://raw.githubusercontent.com/KLemonLab/DpiMGE_Manuscript/master/analysis_MGEs/BOFFO/Selected_Intron/GC_00000001/html/GC_00000001_FWD.html?token=ABXYP4EVFJWYFE2WC6AZZ63ANNNXO).
-   **Real Transposases:** GC\_00000003, GC\_00000040, GC\_00000055,
    GC\_00001693, GC\_00002092, GC\_00002210, GC\_00002310 and
    GC\_00002501.
-   **Integrases (rve domain):** GC\_00000028, GC\_00000085,
    GC\_00001701, GC\_00001775 and GC\_00002348.

We read the BOFFO outputs:

``` r
Intron <- read_excel("analysis_MGEs/Intron.xlsx")
Real_Transposases <- read_excel("analysis_MGEs/Real_Transposases.xlsx")
Integrases_rve <- read_excel("analysis_MGEs/Integrases_rve.xlsx")
```

GCs renamed with MGE type

``` r
Intron <- Intron %>% mutate(MGE = gsub("GC", "INTRON_GC", gene_name))
Real_Transposases <- Real_Transposases %>% mutate(MGE = gsub("GC", "TRANSPOSASE_GC", gene_name))
Integrases_rve <- Integrases_rve %>% mutate(MGE = gsub("GC", "INTEGRASE_GC", gene_name))

ALL <- rbind(Intron, Real_Transposases, Integrases_rve)
```

Summary matrix table with a Genome in each row and each GC listed in
each column. Genomes renamed with full IDs.

``` r
MatrixBOFFO_ALL <- ALL %>%
  group_by(genome_name, MGE) %>%
  summarize(n = n())
MatrixBOFFO_ALL <- spread(MatrixBOFFO_ALL, MGE, n)
MatrixBOFFO_ALL$genomes <- c('ATCC_51524','KPL1914','KPL1922_CDC39_95','KPL1930_CDC2949_98','KPL1931_CDC4294_98','KPL1932_CDC4420_98','KPL1933_CDC4545_98','KPL1934_CDC4709_98','KPL1937_CDC4199_99','KPL1938_CDC4791_99','KPL1939_CDC4792_99','KPL3033','KPL3043','KPL3050','KPL3052','KPL3065','KPL3069','KPL3070','KPL3077','KPL3084','KPL3086','KPL3090','KPL3246','KPL3250','KPL3256','KPL3264','KPL3274','KPL3911')
MatrixBOFFO_ALL <- MatrixBOFFO_ALL %>% remove_rownames %>% column_to_rownames(var="genomes")
MatrixBOFFO_ALL <- select(MatrixBOFFO_ALL, -genome_name)
```

Totals/Stats by column/row:

``` r
GCTotal <- colSums(MatrixBOFFO_ALL, na.rm=TRUE)
GCMean <- round(colMeans(MatrixBOFFO_ALL, na.rm=TRUE), 2)
GCMedian <- round(colMedians(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2)
GCVariance <- round(colVars(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) #Sample variance
GCSD <- round(colSds(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) #Sample standard deviation
GCMAD <- round(colMads(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) #Median absolute deviation
GCMin <- round(colMins(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) 
GCMax <- round(colMaxs(as.matrix(MatrixBOFFO_ALL), na.rm=TRUE), 2) 

MatrixBOFFO_ALL["Total" ,] <- GCTotal
MatrixBOFFO_ALL["Mean" ,] <- GCMean
MatrixBOFFO_ALL["Median" ,] <- GCMedian
MatrixBOFFO_ALL["Variance" ,] <- GCVariance
MatrixBOFFO_ALL["SD" ,] <- GCSD
MatrixBOFFO_ALL["MAD" ,] <- GCMAD
MatrixBOFFO_ALL["Min" ,] <- GCMin
MatrixBOFFO_ALL["Max" ,] <- GCMax
```

Individual tables by MGE:

``` r
Matrix_Intron <- MatrixBOFFO_ALL %>% select(starts_with("INTRON_GC"))
Matrix_Transposases <- MatrixBOFFO_ALL %>% select(starts_with("TRANSPOSASE_GC"))
Matrix_Integrases <- MatrixBOFFO_ALL %>% select(starts_with("INTEGRASE_GC"))

Matrix_Transposases$Total <- rowSums(Matrix_Transposases, na.rm=TRUE)
Matrix_Integrases$Total <- rowSums(Matrix_Integrases, na.rm=TRUE)
```

#### Intron Stats

``` r
opts <- options(knitr.kable.NA = "")
kable(Matrix_Intron, align = "c")
```

|                      | INTRON\_GC\_00000001 |
|:---------------------|:--------------------:|
| ATCC\_51524          |         1.00         |
| KPL1914              |         2.00         |
| KPL1922\_CDC39\_95   |         1.00         |
| KPL1930\_CDC2949\_98 |         2.00         |
| KPL1931\_CDC4294\_98 |         2.00         |
| KPL1932\_CDC4420\_98 |         3.00         |
| KPL1933\_CDC4545\_98 |         1.00         |
| KPL1934\_CDC4709\_98 |        11.00         |
| KPL1937\_CDC4199\_99 |         2.00         |
| KPL1938\_CDC4791\_99 |         1.00         |
| KPL1939\_CDC4792\_99 |         2.00         |
| KPL3033              |         5.00         |
| KPL3043              |         9.00         |
| KPL3050              |        10.00         |
| KPL3052              |         4.00         |
| KPL3065              |         7.00         |
| KPL3069              |         2.00         |
| KPL3070              |         5.00         |
| KPL3077              |         2.00         |
| KPL3084              |         6.00         |
| KPL3086              |         6.00         |
| KPL3090              |        14.00         |
| KPL3246              |         8.00         |
| KPL3250              |         9.00         |
| KPL3256              |         3.00         |
| KPL3264              |         5.00         |
| KPL3274              |         1.00         |
| KPL3911              |         8.00         |
| Total                |        132.00        |
| Mean                 |         4.71         |
| Median               |         3.50         |
| Variance             |        12.88         |
| SD                   |         3.59         |
| MAD                  |         3.71         |
| Min                  |         1.00         |
| Max                  |        14.00         |

#### Transposases Stats

``` r
opts <- options(knitr.kable.NA = "")
kable(Matrix_Transposases, align = "c")
```

|                      | TRANSPOSASE\_GC\_00000003 | TRANSPOSASE\_GC\_00000040 | TRANSPOSASE\_GC\_00000055 | TRANSPOSASE\_GC\_00001693 | TRANSPOSASE\_GC\_00002092 | TRANSPOSASE\_GC\_00002210 | TRANSPOSASE\_GC\_00002310 | TRANSPOSASE\_GC\_00002501 | Total  |
|:---------------------|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:------:|
| ATCC\_51524          |           3.00            |           1.00            |           1.00            |             1             |                           |                           |                           |                           |  6.00  |
| KPL1914              |                           |                           |                           |                           |                           |                           |                           |                           |  0.00  |
| KPL1922\_CDC39\_95   |           1.00            |                           |                           |                           |                           |                           |                           |                           |  1.00  |
| KPL1930\_CDC2949\_98 |                           |                           |                           |                           |                           |                           |                           |                           |  0.00  |
| KPL1931\_CDC4294\_98 |                           |                           |                           |                           |                           |                           |                           |             1             |  1.00  |
| KPL1932\_CDC4420\_98 |           2.00            |                           |                           |                           |                           |                           |                           |                           |  2.00  |
| KPL1933\_CDC4545\_98 |                           |           1.00            |                           |             1             |             1             |             1             |                           |                           |  4.00  |
| KPL1934\_CDC4709\_98 |           6.00            |                           |                           |             1             |                           |                           |                           |                           |  7.00  |
| KPL1937\_CDC4199\_99 |           2.00            |           1.00            |           1.00            |                           |                           |                           |                           |                           |  4.00  |
| KPL1938\_CDC4791\_99 |           1.00            |                           |                           |             1             |                           |                           |                           |                           |  2.00  |
| KPL1939\_CDC4792\_99 |                           |           1.00            |                           |                           |             1             |             1             |                           |                           |  3.00  |
| KPL3033              |           1.00            |                           |                           |                           |                           |                           |                           |                           |  1.00  |
| KPL3043              |           1.00            |                           |                           |                           |                           |                           |                           |                           |  1.00  |
| KPL3050              |           4.00            |                           |                           |             1             |                           |                           |                           |                           |  5.00  |
| KPL3052              |           1.00            |                           |                           |                           |                           |                           |             1             |                           |  2.00  |
| KPL3065              |           1.00            |                           |                           |                           |                           |                           |                           |                           |  1.00  |
| KPL3069              |           1.00            |                           |                           |                           |                           |                           |                           |                           |  1.00  |
| KPL3070              |           5.00            |           1.00            |           1.00            |             1             |                           |                           |                           |                           |  8.00  |
| KPL3077              |           11.00           |           1.00            |           1.00            |                           |                           |                           |                           |                           | 13.00  |
| KPL3084              |           7.00            |           2.00            |           2.00            |             1             |                           |                           |                           |                           | 12.00  |
| KPL3086              |           1.00            |                           |                           |                           |                           |                           |                           |                           |  1.00  |
| KPL3090              |           6.00            |                           |                           |             1             |                           |                           |                           |                           |  7.00  |
| KPL3246              |           3.00            |                           |                           |                           |                           |                           |                           |                           |  3.00  |
| KPL3250              |           2.00            |           1.00            |                           |                           |             1             |             1             |                           |                           |  5.00  |
| KPL3256              |           1.00            |           1.00            |           1.00            |             1             |                           |                           |                           |                           |  4.00  |
| KPL3264              |           6.00            |                           |                           |             1             |                           |                           |                           |                           |  7.00  |
| KPL3274              |                           |                           |                           |                           |                           |                           |             1             |                           |  1.00  |
| KPL3911              |           8.00            |           2.00            |           2.00            |             1             |                           |                           |                           |                           | 13.00  |
| Total                |           74.00           |           12.00           |           9.00            |            11             |             3             |             3             |             2             |             1             | 115.00 |
| Mean                 |           3.36            |           1.20            |           1.29            |             1             |             1             |             1             |             1             |             1             | 10.85  |
| Median               |           2.00            |           1.00            |           1.00            |             1             |             1             |             1             |             1             |             1             |  9.00  |
| Variance             |           8.24            |           0.18            |           0.24            |             0             |             0             |             0             |             0             |                           |  8.66  |
| SD                   |           2.87            |           0.42            |           0.49            |             0             |             0             |             0             |             0             |                           |  3.78  |
| MAD                  |           1.48            |           0.00            |           0.00            |             0             |             0             |             0             |             0             |             0             |  1.48  |
| Min                  |           1.00            |           1.00            |           1.00            |             1             |             1             |             1             |             1             |             1             |  8.00  |
| Max                  |           11.00           |           2.00            |           2.00            |             1             |             1             |             1             |             1             |             1             | 20.00  |

#### Integrases Stats

``` r
opts <- options(knitr.kable.NA = "")
kable(Matrix_Integrases, align = "c")
```

|                      | INTEGRASE\_GC\_00000028 | INTEGRASE\_GC\_00000085 | INTEGRASE\_GC\_00001701 | INTEGRASE\_GC\_00001775 | INTEGRASE\_GC\_00002348 | Total |
|:---------------------|:-----------------------:|:-----------------------:|:-----------------------:|:-----------------------:|:-----------------------:|:-----:|
| ATCC\_51524          |                         |            1            |                         |                         |                         | 1.00  |
| KPL1914              |                         |            1            |                         |            1            |                         | 2.00  |
| KPL1922\_CDC39\_95   |                         |                         |                         |                         |                         | 0.00  |
| KPL1930\_CDC2949\_98 |          1.00           |                         |                         |                         |                         | 1.00  |
| KPL1931\_CDC4294\_98 |                         |                         |                         |                         |                         | 0.00  |
| KPL1932\_CDC4420\_98 |                         |                         |          1.00           |            1            |                         | 2.00  |
| KPL1933\_CDC4545\_98 |          2.00           |                         |                         |                         |                         | 2.00  |
| KPL1934\_CDC4709\_98 |                         |            1            |          1.00           |                         |                         | 2.00  |
| KPL1937\_CDC4199\_99 |                         |                         |                         |                         |                         | 0.00  |
| KPL1938\_CDC4791\_99 |                         |            1            |                         |                         |                         | 1.00  |
| KPL1939\_CDC4792\_99 |          1.00           |                         |                         |                         |                         | 1.00  |
| KPL3033              |                         |            1            |                         |                         |                         | 1.00  |
| KPL3043              |          1.00           |                         |                         |                         |                         | 1.00  |
| KPL3050              |          2.00           |            1            |          3.00           |                         |                         | 6.00  |
| KPL3052              |                         |            1            |          1.00           |                         |            1            | 3.00  |
| KPL3065              |          1.00           |                         |                         |                         |                         | 1.00  |
| KPL3069              |          1.00           |                         |          1.00           |                         |                         | 2.00  |
| KPL3070              |                         |            1            |                         |                         |                         | 1.00  |
| KPL3077              |          1.00           |                         |                         |            1            |                         | 2.00  |
| KPL3084              |                         |            1            |                         |                         |                         | 1.00  |
| KPL3086              |                         |                         |                         |                         |                         | 0.00  |
| KPL3090              |          1.00           |            1            |          1.00           |                         |                         | 3.00  |
| KPL3246              |                         |                         |                         |                         |                         | 0.00  |
| KPL3250              |          1.00           |            1            |                         |                         |                         | 2.00  |
| KPL3256              |                         |                         |                         |                         |                         | 0.00  |
| KPL3264              |                         |            1            |                         |                         |                         | 1.00  |
| KPL3274              |          5.00           |                         |                         |                         |            1            | 6.00  |
| KPL3911              |                         |            1            |                         |                         |                         | 1.00  |
| Total                |          17.00          |           13            |          8.00           |            3            |            2            | 43.00 |
| Mean                 |          1.55           |            1            |          1.33           |            1            |            1            | 5.88  |
| Median               |          1.00           |            1            |          1.00           |            1            |            1            | 5.00  |
| Variance             |          1.47           |            0            |          0.67           |            0            |            0            | 2.14  |
| SD                   |          1.21           |            0            |          0.82           |            0            |            0            | 2.03  |
| MAD                  |          0.00           |            0            |          0.00           |            0            |            0            | 0.00  |
| Min                  |          1.00           |            1            |          1.00           |            1            |            1            | 5.00  |
| Max                  |          5.00           |            1            |          3.00           |            1            |            1            | 11.00 |

``` r
write.csv(Matrix_Intron, "analysis_MGEs/Intron_Table.csv", row.names = TRUE)
write.csv(Matrix_Transposases, "analysis_MGEs/Transposases_Table.csv", row.names = TRUE)
write.csv(Matrix_Integrases, "analysis_MGEs/Integrases_Table.csv", row.names = TRUE)
```

# <u>REFERENCES</u>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Galperin2021" class="csl-entry">

Galperin, Michael Y, Yuri I Wolf, Kira S Makarova, Roberto Vera Alvarez,
David Landsman, and Eugene V Koonin. 2020. “COG Database Update: Focus
on Microbial Diversity, Model Organisms, and Widespread Pathogens.”
*Nucleic Acids Research* 49 (D1): D274–81.
<https://doi.org/10.1093/nar/gkaa1018>.

</div>

<div id="ref-Gilchrist2021" class="csl-entry">

Gilchrist, Rachel J., Lisa M. Gunter, Samantha F. Anderson, and Clive
D.L. Wynne. 2021. “The Click Is Not the Trick: The Efficacy of Clickers
and Other Reinforcement Methods in Training Naïve Dogs to Perform New
Tasks.” *PeerJ* 9 (February): e10881.
<https://doi.org/10.7717/peerj.10881>.

</div>

<div id="ref-kanehisa2000" class="csl-entry">

Kanehisa, M. 2000. “KEGG: Kyoto Encyclopedia of Genes and Genomes.”
*Nucleic Acids Research* 28 (1): 27–30.
<https://doi.org/10.1093/nar/28.1.27>.

</div>

<div id="ref-kanehisa2015" class="csl-entry">

Kanehisa, Minoru, Yoko Sato, Masayuki Kawashima, Miho Furumichi, and Mao
Tanabe. 2015. “KEGG as a Reference Resource for Gene and Protein
Annotation.” *Nucleic Acids Research* 44 (D1): D457–62.
<https://doi.org/10.1093/nar/gkv1070>.

</div>

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

<div id="ref-Seemann2014" class="csl-entry">

Seemann, T. 2014. “Prokka: Rapid Prokaryotic Genome Annotation.”
*Bioinformatics* 30 (14): 2068–69.
<https://doi.org/10.1093/bioinformatics/btu153>.

</div>

<div id="ref-Tatusov1997" class="csl-entry">

Tatusov, R. L. 1997. “A Genomic Perspective on Protein Families.”
*Science* 278 (5338): 631–37.
<https://doi.org/10.1126/science.278.5338.631>.

</div>

</div>
