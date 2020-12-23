*Dolosigranulum pigrum* COG ANALYSIS
================

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

\*“Core” is used in the code to avoid problems with the “/” symbol and
later replaced with “Soft/Core” for plotting.

There are 1515 gene clusters (GC) (52.2%) in the Accessory vs. 1388
(47.8%) in the Soft/Core at the pangenome level

## COG Analysis at the Gene Level

We define a new variable `COGs` to use in the plots. This variable is
based on `COG_CATEGORY` but with a cleaner definition of unclassified,
uninformative or mixed assignments: \* COG categories “Function Unknown”
and “General function predictions only” were considered as
“Uninformative”. \* If the COG category is a mix (for example G|S|M)
it gets labeled like “Ambiguous”. \* Also missing values (NA) are
labeled as Unclassified".

``` r
DpigPangenome$COGs <- DpigPangenome$COG_CATEGORY
DpigPangenome$COGs[DpigPangenome$COGs =="S"]<- "Uninformative"
DpigPangenome$COGs[DpigPangenome$COGs =="R"]<- "Uninformative"
DpigPangenome$COGs[grepl('|', DpigPangenome$COGs,fixed=TRUE)]<-"Ambiguous"
DpigPangenome$COGs[is.na(DpigPangenome$COGs)]<-"Unclassified"
```

Summary of GOC annotated genes:

<table>

<thead>

<tr>

<th style="text-align:left;">

Genes

</th>

<th style="text-align:right;">

Count

</th>

<th style="text-align:right;">

Percentage

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Total in Dpig Pangenome

</td>

<td style="text-align:right;">

49415

</td>

<td style="text-align:right;">

100.0

</td>

</tr>

<tr>

<td style="text-align:left;">

COG Category Uninformative = Function Unknown

</td>

<td style="text-align:right;">

2745

</td>

<td style="text-align:right;">

5.6

</td>

</tr>

<tr>

<td style="text-align:left;">

COG Category Uninformative = General function predictions only

</td>

<td style="text-align:right;">

2833

</td>

<td style="text-align:right;">

5.7

</td>

</tr>

<tr>

<td style="text-align:left;">

COG Category Ambiguous (Mixed COG Category)

</td>

<td style="text-align:right;">

4953

</td>

<td style="text-align:right;">

10.0

</td>

</tr>

<tr>

<td style="text-align:left;">

COG Category Unclassified (Non-assigned)

</td>

<td style="text-align:right;">

8378

</td>

<td style="text-align:right;">

17.0

</td>

</tr>

<tr>

<td style="text-align:left;">

Informative COGs (Total - Uninformative, Ambiguous & Unclassified)

</td>

<td style="text-align:right;">

30506

</td>

<td style="text-align:right;">

61.7

</td>

</tr>

</tbody>

</table>

60.7% of the gene calls are Informative.

## COG Analysis at the Gene Cluster Level

This analysis was done at the pangenomic gene cluster level (GC). Since
many gene clusters had mixed COG category assignments a solution is to
assign each individual gene call to their corresponding
Genome/AccessoryvsCore/COG grouping weighting their contribution by
dividing their count by the number of genes in their GC.

### GCs by COG Category and Genome

The table “GCsbyCOG\_Genome” groups the genes by genome; and inside
genomes by Accessory vs. Soft/Core status, and nested inside as the COG
category. But, in this case, instead of counting the elements in each
group we calculated the sum of 1/num\_genes\_in\_gene\_cluster.

``` r
GCsbyCOG_Genome <- DpigPangenome %>%
  group_by(genome_name, accessory_vs_core, COGs) %>%
  summarise(num_corrected_genes=sum(1/num_genes_in_gene_cluster))
```

The total sum of all values in the `num_corrected_genes` variable should
add up to the number of CGs:

``` r
sum(GCsbyCOG_Genome$num_corrected_genes)
```

    ## [1] 2903

``` r
nrow(DpigPangenome %>% group_by(gene_cluster_id) %>% summarise)
```

    ## [1] 2903

Extra column to label the gray scale portion of the plots:

``` r
GCsbyCOG_Genome <- GCsbyCOG_Genome %>%
  mutate(Assignment=ifelse(COGs!="Uninformative" & COGs!="Ambiguous" & COGs!="Unclassified", "Informative", COGs))
```

Summary of GOC annotated GCs in the Accessory vs. Soft/Core :

``` r
TableGC <- GCsbyCOG_Genome %>% 
  group_by(accessory_vs_core, Assignment) %>%
  summarize(corrected_genes=sum(num_corrected_genes))

TableGC$Percentages <- round(100*TableGC$corrected_genes/sum(TableGC$corrected_genes), 1)

kable(TableGC)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

accessory\_vs\_core

</th>

<th style="text-align:left;">

Assignment

</th>

<th style="text-align:right;">

corrected\_genes

</th>

<th style="text-align:right;">

Percentages

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Ambiguous

</td>

<td style="text-align:right;">

105.09970

</td>

<td style="text-align:right;">

3.6

</td>

</tr>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Informative

</td>

<td style="text-align:right;">

540.43075

</td>

<td style="text-align:right;">

18.6

</td>

</tr>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Unclassified

</td>

<td style="text-align:right;">

809.50850

</td>

<td style="text-align:right;">

27.9

</td>

</tr>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Uninformative

</td>

<td style="text-align:right;">

59.96105

</td>

<td style="text-align:right;">

2.1

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Ambiguous

</td>

<td style="text-align:right;">

142.29524

</td>

<td style="text-align:right;">

4.9

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Informative

</td>

<td style="text-align:right;">

922.49256

</td>

<td style="text-align:right;">

31.8

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Unclassified

</td>

<td style="text-align:right;">

141.07739

</td>

<td style="text-align:right;">

4.9

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Uninformative

</td>

<td style="text-align:right;">

182.13481

</td>

<td style="text-align:right;">

6.3

</td>

</tr>

</tbody>

</table>

Summary of GOC annotated GCs in the Accessory:

``` r
TableGCAccessory <- GCsbyCOG_Genome %>% 
  filter(accessory_vs_core =="Accessory") %>%
  group_by(accessory_vs_core, Assignment) %>%
  summarize(corrected_genes=sum(num_corrected_genes))

TableGCAccessory$Percentages <- round(100*TableGCAccessory$corrected_genes/sum(TableGCAccessory$corrected_genes), 1)

kable(TableGCAccessory)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

accessory\_vs\_core

</th>

<th style="text-align:left;">

Assignment

</th>

<th style="text-align:right;">

corrected\_genes

</th>

<th style="text-align:right;">

Percentages

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Ambiguous

</td>

<td style="text-align:right;">

105.09970

</td>

<td style="text-align:right;">

6.9

</td>

</tr>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Informative

</td>

<td style="text-align:right;">

540.43075

</td>

<td style="text-align:right;">

35.7

</td>

</tr>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Unclassified

</td>

<td style="text-align:right;">

809.50850

</td>

<td style="text-align:right;">

53.4

</td>

</tr>

<tr>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:left;">

Uninformative

</td>

<td style="text-align:right;">

59.96105

</td>

<td style="text-align:right;">

4.0

</td>

</tr>

</tbody>

</table>

Summary of GOC annotated GCs in the Soft/Core :

``` r
TableGCCore <- GCsbyCOG_Genome %>% 
  filter(accessory_vs_core =="Core") %>%
  group_by(accessory_vs_core, Assignment) %>%
  summarize(corrected_genes=sum(num_corrected_genes))

TableGCCore$Percentages <- round(100*TableGCCore$corrected_genes/sum(TableGCCore$corrected_genes), 1)

kable(TableGCCore)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

accessory\_vs\_core

</th>

<th style="text-align:left;">

Assignment

</th>

<th style="text-align:right;">

corrected\_genes

</th>

<th style="text-align:right;">

Percentages

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Ambiguous

</td>

<td style="text-align:right;">

142.2952

</td>

<td style="text-align:right;">

10.3

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Informative

</td>

<td style="text-align:right;">

922.4926

</td>

<td style="text-align:right;">

66.5

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Unclassified

</td>

<td style="text-align:right;">

141.0774

</td>

<td style="text-align:right;">

10.2

</td>

</tr>

<tr>

<td style="text-align:left;">

Core

</td>

<td style="text-align:left;">

Uninformative

</td>

<td style="text-align:right;">

182.1348

</td>

<td style="text-align:right;">

13.1

</td>

</tr>

</tbody>

</table>

Summary of GOC annotated GCs by Genome in the Accessory vs. Soft/Core :

``` r
TableGenomes <- GCsbyCOG_Genome %>% 
  group_by(genome_name, accessory_vs_core) %>% 
  summarize(corrected_genes=sum(num_corrected_genes))

kable(TableGenomes)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

genome\_name

</th>

<th style="text-align:left;">

accessory\_vs\_core

</th>

<th style="text-align:right;">

corrected\_genes

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

ATCC\_51524

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

46.24398

</td>

</tr>

<tr>

<td style="text-align:left;">

ATCC\_51524

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.26676

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1914

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

44.52563

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1914

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.58011

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1922\_CDC39\_95

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

59.17963

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1922\_CDC39\_95

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.27983

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1930\_CDC2949\_98

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

48.40225

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1930\_CDC2949\_98

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.09844

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1931\_CDC4294\_98

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

180.28196

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1931\_CDC4294\_98

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.22178

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1932\_CDC4420\_98

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

70.59157

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1932\_CDC4420\_98

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.44550

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1933\_CDC4545\_98

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

53.11917

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1933\_CDC4545\_98

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.18986

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1934\_CDC4709\_98

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

35.73283

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1934\_CDC4709\_98

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

50.82800

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1937\_CDC4199\_99

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

72.05543

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1937\_CDC4199\_99

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

48.69897

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1938\_CDC4791\_99

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

50.05558

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1938\_CDC4791\_99

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.15724

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1939\_CDC4792\_99

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

65.32743

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL1939\_CDC4792\_99

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.42181

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3033

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

69.20019

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3033

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.42137

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3043

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

27.79795

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3043

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.53672

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3050

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

58.42861

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3050

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.80908

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3052

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

49.69199

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3052

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.84545

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3065

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

22.35714

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3065

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.67527

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3069

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

77.61282

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3069

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

50.02098

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3070

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

28.60126

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3070

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.94724

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3077

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

55.70137

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3077

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.61444

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3084

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

30.26084

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3084

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.81861

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3086

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

24.66277

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3086

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.87085

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3090

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

86.80387

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3090

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.81137

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3246

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

32.51061

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3246

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.63258

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3250

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

33.54099

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3250

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.59168

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3256

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

68.10628

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3256

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.35919

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3264

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

58.50242

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3264

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.70653

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3274

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

35.04517

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3274

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.18935

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3911

</td>

<td style="text-align:left;">

Accessory

</td>

<td style="text-align:right;">

30.66027

</td>

</tr>

<tr>

<td style="text-align:left;">

KPL3911

</td>

<td style="text-align:left;">

Core

</td>

<td style="text-align:right;">

49.96100

</td>

</tr>

</tbody>

</table>

Renaming and ordering variables factor levels for plotting:

``` r
GCsbyCOG_Genome$accessory_vs_core <- factor(GCsbyCOG_Genome$accessory_vs_core, levels =c("Core", "Accessory"))

GCsbyCOG_Genome$COGs <- recode_factor(GCsbyCOG_Genome$COGs, "Q"="Secondary metabolites biosynthesis, transport, and catabolism","P"="Inorganic ion transport and metabolism","I"="Lipid transport and metabolism","H"="Coenzyme transport and metabolism","G"="Carbohydrate transport and metabolism","F"="Nucleotide transport and metabolism","E"="Amino acid transport and metabolism","C"="Energy production and conversion","X"="Mobilome: prophages, transposons","L"="Replication, recombination and repair","K"="Transcription","J"="Translation, ribosomal structure and biogenesis","V"="Defense mechanisms","U"="Intracellular trafficking, secretion, and vesicular transport","T"="Signal transduction mechanisms","O"="Post-translational modification, protein turnover, and chaperones","N"="Cell Motility","M"="Cell wall/membrane/envelope biogenesis","D"="Cell cycle control, cell division, chromosome partitioning","Uninformative"="Uninformative","Ambiguous"="Ambiguous","Unclassified"="Unclassified", .ordered = TRUE)

GCsbyCOG_Genome$Assignment <- recode_factor(GCsbyCOG_Genome$Assignment,  "Informative"=" ", "Uninformative"="Uninformative", "Ambiguous"="Ambiguous", "Unclassified"="Unclassified", .ordered = TRUE)

GCsbyCOG_Genome$genome_name <- recode_factor(GCsbyCOG_Genome$genome_name, "Dpigrum_ATCC_51524"="ATCC 51524", "Dpigrum_KPL3250"="KPL3250", "Dpigrum_KPL1939_CDC4792_99"="CDC 4792-99","Dpigrum_KPL1934_CDC4709_98"="CDC 4709-98", "Dpigrum_KPL1922_CDC39_95"="CDC 39-95", "Dpigrum_KPL3264"="KPL3264", "Dpigrum_KPL3256"="KPL3256", "Dpigrum_KPL3033"="KPL3033", "Dpigrum_KPL1933_CDC4545_98"="CDC 4545-98", "Dpigrum_KPL1930_CDC2949_98"="CDC 2949-98", "Dpigrum_KPL3069"="KPL3069", "Dpigrum_KPL3052"="KPL3052", "Dpigrum_KPL3090"="KPL3090", "Dpigrum_KPL3086"="KPL3086", "Dpigrum_KPL3065"="KPL3065", "Dpigrum_KPL3043"="KPL3043", "Dpigrum_KPL3911"="KPL3911", "Dpigrum_KPL3084"="KPL3084", "Dpigrum_KPL3070"="KPL3070",
"Dpigrum_KPL3246"="KPL3246", "Dpigrum_KPL1937_CDC4199_99"="CDC 4199-99","Dpigrum_KPL3274"="KPL3274","Dpigrum_KPL3050"="KPL3050","Dpigrum_KPL1938_CDC4791_99"="CDC 4791-99", "Dpigrum_KPL1932_CDC4420_98"="CDC 4420-98", "Dpigrum_KPL3077"="KPL3077", "Dpigrum_KPL1931_CDC4294_98"="CDC 4294-98", "Dpigrum_KPL1914"="KPL1914", .ordered = TRUE)
```

### GCs by COG Category

The table “GCsbyCOG” groups the genes by Accessory vs. Soft/Core status,
and nested inside as the COG category.

``` r
GCsbyCOG <- DpigPangenome %>%
  group_by(accessory_vs_core, COGs) %>%
  summarise(num_corrected_genes=sum(1/num_genes_in_gene_cluster))
```

Renaming and ordering variables factor levels for plotting:

``` r
GCsbyCOG$COGs <- recode_factor(GCsbyCOG$COGs, "Q"="Secondary metabolites biosynthesis, transport, and catabolism","P"="Inorganic ion transport and metabolism","I"="Lipid transport and metabolism","H"="Coenzyme transport and metabolism","G"="Carbohydrate transport and metabolism","F"="Nucleotide transport and metabolism","E"="Amino acid transport and metabolism","C"="Energy production and conversion","X"="Mobilome: prophages, transposons","L"="Replication, recombination and repair","K"="Transcription","J"="Translation, ribosomal structure and biogenesis","V"="Defense mechanisms","U"="Intracellular trafficking, secretion, and vesicular transport","T"="Signal transduction mechanisms","O"="Post-translational modification, protein turnover, and chaperones","N"="Cell Motility","M"="Cell wall/membrane/envelope biogenesis","D"="Cell cycle control, cell division, chromosome partitioning","Uninformative"="Uninformative","Ambiguous"="Ambiguous","Unclassified"="Unclassified", .ordered = TRUE)
```

New table “GCsbyCOG\_CorevsAcc” in wide format. % of each category
relative to the “Accessory” or “Soft/Core” calculated (pTotal.
variables). Also, total GCs for each COG category calculated, as well as
% of GCs in the “Accessory” and “Soft/Core” relative to each category
(p. values). The ratio between the number of GC in the “Accessory”
vs. the “Soft/Core” is calculated for each COG:

``` r
GCsbyCOG_CorevsAcc <- spread(GCsbyCOG, accessory_vs_core, num_corrected_genes)
GCsbyCOG_CorevsAcc$pTotal.Accessory <- round(100*GCsbyCOG_CorevsAcc$Accessory/sum(GCsbyCOG_CorevsAcc$Accessory), 1)
GCsbyCOG_CorevsAcc$pTotal.Core <- round(100*GCsbyCOG_CorevsAcc$Core/sum(GCsbyCOG_CorevsAcc$Core), 1)
GCsbyCOG_CorevsAcc$total <- GCsbyCOG_CorevsAcc$Accessory + GCsbyCOG_CorevsAcc$Core
GCsbyCOG_CorevsAcc$pTotal.total <- round(100*GCsbyCOG_CorevsAcc$total/sum(GCsbyCOG_CorevsAcc$total), 1)
GCsbyCOG_CorevsAcc$p.accessory <- round(100*(GCsbyCOG_CorevsAcc$Accessory/GCsbyCOG_CorevsAcc$total), 1)
GCsbyCOG_CorevsAcc$p.core <- round(100*(GCsbyCOG_CorevsAcc$Core/GCsbyCOG_CorevsAcc$total), 1)
GCsbyCOG_CorevsAcc$ratio <- round(GCsbyCOG_CorevsAcc$Accessory/GCsbyCOG_CorevsAcc$Core, 2)

kable(GCsbyCOG_CorevsAcc)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

COGs

</th>

<th style="text-align:right;">

Accessory

</th>

<th style="text-align:right;">

Core

</th>

<th style="text-align:right;">

pTotal.Accessory

</th>

<th style="text-align:right;">

pTotal.Core

</th>

<th style="text-align:right;">

total

</th>

<th style="text-align:right;">

pTotal.total

</th>

<th style="text-align:right;">

p.accessory

</th>

<th style="text-align:right;">

p.core

</th>

<th style="text-align:right;">

ratio

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Secondary metabolites biosynthesis, transport, and catabolism

</td>

<td style="text-align:right;">

1.000000

</td>

<td style="text-align:right;">

5.000000

</td>

<td style="text-align:right;">

0.1

</td>

<td style="text-align:right;">

0.4

</td>

<td style="text-align:right;">

6.000000

</td>

<td style="text-align:right;">

0.2

</td>

<td style="text-align:right;">

16.7

</td>

<td style="text-align:right;">

83.3

</td>

<td style="text-align:right;">

0.20

</td>

</tr>

<tr>

<td style="text-align:left;">

Inorganic ion transport and metabolism

</td>

<td style="text-align:right;">

11.800000

</td>

<td style="text-align:right;">

65.609368

</td>

<td style="text-align:right;">

0.8

</td>

<td style="text-align:right;">

4.7

</td>

<td style="text-align:right;">

77.409368

</td>

<td style="text-align:right;">

2.7

</td>

<td style="text-align:right;">

15.2

</td>

<td style="text-align:right;">

84.8

</td>

<td style="text-align:right;">

0.18

</td>

</tr>

<tr>

<td style="text-align:left;">

Lipid transport and metabolism

</td>

<td style="text-align:right;">

3.000000

</td>

<td style="text-align:right;">

36.866725

</td>

<td style="text-align:right;">

0.2

</td>

<td style="text-align:right;">

2.7

</td>

<td style="text-align:right;">

39.866725

</td>

<td style="text-align:right;">

1.4

</td>

<td style="text-align:right;">

7.5

</td>

<td style="text-align:right;">

92.5

</td>

<td style="text-align:right;">

0.08

</td>

</tr>

<tr>

<td style="text-align:left;">

Coenzyme transport and metabolism

</td>

<td style="text-align:right;">

11.947368

</td>

<td style="text-align:right;">

45.034483

</td>

<td style="text-align:right;">

0.8

</td>

<td style="text-align:right;">

3.2

</td>

<td style="text-align:right;">

56.981851

</td>

<td style="text-align:right;">

2.0

</td>

<td style="text-align:right;">

21.0

</td>

<td style="text-align:right;">

79.0

</td>

<td style="text-align:right;">

0.27

</td>

</tr>

<tr>

<td style="text-align:left;">

Carbohydrate transport and metabolism

</td>

<td style="text-align:right;">

135.146781

</td>

<td style="text-align:right;">

80.606549

</td>

<td style="text-align:right;">

8.9

</td>

<td style="text-align:right;">

5.8

</td>

<td style="text-align:right;">

215.753331

</td>

<td style="text-align:right;">

7.4

</td>

<td style="text-align:right;">

62.6

</td>

<td style="text-align:right;">

37.4

</td>

<td style="text-align:right;">

1.68

</td>

</tr>

<tr>

<td style="text-align:left;">

Nucleotide transport and metabolism

</td>

<td style="text-align:right;">

6.000000

</td>

<td style="text-align:right;">

48.000000

</td>

<td style="text-align:right;">

0.4

</td>

<td style="text-align:right;">

3.5

</td>

<td style="text-align:right;">

54.000000

</td>

<td style="text-align:right;">

1.9

</td>

<td style="text-align:right;">

11.1

</td>

<td style="text-align:right;">

88.9

</td>

<td style="text-align:right;">

0.12

</td>

</tr>

<tr>

<td style="text-align:left;">

Amino acid transport and metabolism

</td>

<td style="text-align:right;">

21.601708

</td>

<td style="text-align:right;">

72.948276

</td>

<td style="text-align:right;">

1.4

</td>

<td style="text-align:right;">

5.3

</td>

<td style="text-align:right;">

94.549984

</td>

<td style="text-align:right;">

3.3

</td>

<td style="text-align:right;">

22.8

</td>

<td style="text-align:right;">

77.2

</td>

<td style="text-align:right;">

0.30

</td>

</tr>

<tr>

<td style="text-align:left;">

Energy production and conversion

</td>

<td style="text-align:right;">

10.558824

</td>

<td style="text-align:right;">

46.536946

</td>

<td style="text-align:right;">

0.7

</td>

<td style="text-align:right;">

3.4

</td>

<td style="text-align:right;">

57.095769

</td>

<td style="text-align:right;">

2.0

</td>

<td style="text-align:right;">

18.5

</td>

<td style="text-align:right;">

81.5

</td>

<td style="text-align:right;">

0.23

</td>

</tr>

<tr>

<td style="text-align:left;">

Mobilome: prophages, transposons

</td>

<td style="text-align:right;">

54.236979

</td>

<td style="text-align:right;">

4.000000

</td>

<td style="text-align:right;">

3.6

</td>

<td style="text-align:right;">

0.3

</td>

<td style="text-align:right;">

58.236979

</td>

<td style="text-align:right;">

2.0

</td>

<td style="text-align:right;">

93.1

</td>

<td style="text-align:right;">

6.9

</td>

<td style="text-align:right;">

13.56

</td>

</tr>

<tr>

<td style="text-align:left;">

Replication, recombination and repair

</td>

<td style="text-align:right;">

58.810201

</td>

<td style="text-align:right;">

72.476833

</td>

<td style="text-align:right;">

3.9

</td>

<td style="text-align:right;">

5.2

</td>

<td style="text-align:right;">

131.287034

</td>

<td style="text-align:right;">

4.5

</td>

<td style="text-align:right;">

44.8

</td>

<td style="text-align:right;">

55.2

</td>

<td style="text-align:right;">

0.81

</td>

</tr>

<tr>

<td style="text-align:left;">

Transcription

</td>

<td style="text-align:right;">

59.731944

</td>

<td style="text-align:right;">

52.312338

</td>

<td style="text-align:right;">

3.9

</td>

<td style="text-align:right;">

3.8

</td>

<td style="text-align:right;">

112.044282

</td>

<td style="text-align:right;">

3.9

</td>

<td style="text-align:right;">

53.3

</td>

<td style="text-align:right;">

46.7

</td>

<td style="text-align:right;">

1.14

</td>

</tr>

<tr>

<td style="text-align:left;">

Translation, ribosomal structure and biogenesis

</td>

<td style="text-align:right;">

18.327500

</td>

<td style="text-align:right;">

171.438894

</td>

<td style="text-align:right;">

1.2

</td>

<td style="text-align:right;">

12.4

</td>

<td style="text-align:right;">

189.766394

</td>

<td style="text-align:right;">

6.5

</td>

<td style="text-align:right;">

9.7

</td>

<td style="text-align:right;">

90.3

</td>

<td style="text-align:right;">

0.11

</td>

</tr>

<tr>

<td style="text-align:left;">

Defense mechanisms

</td>

<td style="text-align:right;">

96.981384

</td>

<td style="text-align:right;">

32.686978

</td>

<td style="text-align:right;">

6.4

</td>

<td style="text-align:right;">

2.4

</td>

<td style="text-align:right;">

129.668361

</td>

<td style="text-align:right;">

4.5

</td>

<td style="text-align:right;">

74.8

</td>

<td style="text-align:right;">

25.2

</td>

<td style="text-align:right;">

2.97

</td>

</tr>

<tr>

<td style="text-align:left;">

Intracellular trafficking, secretion, and vesicular transport

</td>

<td style="text-align:right;">

5.000000

</td>

<td style="text-align:right;">

10.000000

</td>

<td style="text-align:right;">

0.3

</td>

<td style="text-align:right;">

0.7

</td>

<td style="text-align:right;">

15.000000

</td>

<td style="text-align:right;">

0.5

</td>

<td style="text-align:right;">

33.3

</td>

<td style="text-align:right;">

66.7

</td>

<td style="text-align:right;">

0.50

</td>

</tr>

<tr>

<td style="text-align:left;">

Signal transduction mechanisms

</td>

<td style="text-align:right;">

5.351191

</td>

<td style="text-align:right;">

35.696591

</td>

<td style="text-align:right;">

0.4

</td>

<td style="text-align:right;">

2.6

</td>

<td style="text-align:right;">

41.047781

</td>

<td style="text-align:right;">

1.4

</td>

<td style="text-align:right;">

13.0

</td>

<td style="text-align:right;">

87.0

</td>

<td style="text-align:right;">

0.15

</td>

</tr>

<tr>

<td style="text-align:left;">

Post-translational modification, protein turnover, and chaperones

</td>

<td style="text-align:right;">

11.807692

</td>

<td style="text-align:right;">

56.966667

</td>

<td style="text-align:right;">

0.8

</td>

<td style="text-align:right;">

4.1

</td>

<td style="text-align:right;">

68.774359

</td>

<td style="text-align:right;">

2.4

</td>

<td style="text-align:right;">

17.2

</td>

<td style="text-align:right;">

82.8

</td>

<td style="text-align:right;">

0.21

</td>

</tr>

<tr>

<td style="text-align:left;">

Cell Motility

</td>

<td style="text-align:right;">

5.941177

</td>

<td style="text-align:right;">

2.214286

</td>

<td style="text-align:right;">

0.4

</td>

<td style="text-align:right;">

0.2

</td>

<td style="text-align:right;">

8.155462

</td>

<td style="text-align:right;">

0.3

</td>

<td style="text-align:right;">

72.8

</td>

<td style="text-align:right;">

27.2

</td>

<td style="text-align:right;">

2.68

</td>

</tr>

<tr>

<td style="text-align:left;">

Cell wall/membrane/envelope biogenesis

</td>

<td style="text-align:right;">

15.584827

</td>

<td style="text-align:right;">

60.917824

</td>

<td style="text-align:right;">

1.0

</td>

<td style="text-align:right;">

4.4

</td>

<td style="text-align:right;">

76.502651

</td>

<td style="text-align:right;">

2.6

</td>

<td style="text-align:right;">

20.4

</td>

<td style="text-align:right;">

79.6

</td>

<td style="text-align:right;">

0.26

</td>

</tr>

<tr>

<td style="text-align:left;">

Cell cycle control, cell division, chromosome partitioning

</td>

<td style="text-align:right;">

7.603175

</td>

<td style="text-align:right;">

23.179803

</td>

<td style="text-align:right;">

0.5

</td>

<td style="text-align:right;">

1.7

</td>

<td style="text-align:right;">

30.782978

</td>

<td style="text-align:right;">

1.1

</td>

<td style="text-align:right;">

24.7

</td>

<td style="text-align:right;">

75.3

</td>

<td style="text-align:right;">

0.33

</td>

</tr>

<tr>

<td style="text-align:left;">

Uninformative

</td>

<td style="text-align:right;">

59.961053

</td>

<td style="text-align:right;">

182.134812

</td>

<td style="text-align:right;">

4.0

</td>

<td style="text-align:right;">

13.1

</td>

<td style="text-align:right;">

242.095865

</td>

<td style="text-align:right;">

8.3

</td>

<td style="text-align:right;">

24.8

</td>

<td style="text-align:right;">

75.2

</td>

<td style="text-align:right;">

0.33

</td>

</tr>

<tr>

<td style="text-align:left;">

Ambiguous

</td>

<td style="text-align:right;">

105.099695

</td>

<td style="text-align:right;">

142.295243

</td>

<td style="text-align:right;">

6.9

</td>

<td style="text-align:right;">

10.3

</td>

<td style="text-align:right;">

247.394938

</td>

<td style="text-align:right;">

8.5

</td>

<td style="text-align:right;">

42.5

</td>

<td style="text-align:right;">

57.5

</td>

<td style="text-align:right;">

0.74

</td>

</tr>

<tr>

<td style="text-align:left;">

Unclassified

</td>

<td style="text-align:right;">

809.508501

</td>

<td style="text-align:right;">

141.077387

</td>

<td style="text-align:right;">

53.4

</td>

<td style="text-align:right;">

10.2

</td>

<td style="text-align:right;">

950.585889

</td>

<td style="text-align:right;">

32.7

</td>

<td style="text-align:right;">

85.2

</td>

<td style="text-align:right;">

14.8

</td>

<td style="text-align:right;">

5.74

</td>

</tr>

</tbody>

</table>

## Plots

Color Palettes

``` r
getPalette <- colorRampPalette(brewer.pal(8, "Set1"))
CountTotalCOGs <- length(unique(GCsbyCOG_Genome$COGs))

palette1 <- c("grey60", "grey40", "grey20", getPalette(CountTotalCOGs-3)) # 22 elements: Colors + Grays
palette2 <- getPalette(CountTotalCOGs-3) # 19 elements: Colors
palette3 <- c("grey60", "grey40", "grey20", "white") # 4 elements: White + Grays
```

### Plots Accessory vs. Soft/Core

Panel A in main figure:

``` r
pA <- ggplot(GCsbyCOG_Genome, aes(x = accessory_vs_core, y = num_corrected_genes, fill = fct_rev(COGs))) +
  stat_summary(fun=sum ,geom="bar", position = "stack") +
  scale_x_discrete(labels = c("Soft/Core", "Accessory")) +
  scale_fill_manual(values = palette1) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 1500, by = 250)) +
  labs(fill="COG Categories", x=" ", y= "Number of Gene Clusters") +
  theme_classic() +
  theme(axis.title = element_text(size = 9), axis.text = element_text(size=7), plot.margin=unit(c(10,0,10,20),"pt"), legend.position = "none") 
```

This plot is used for the grayscale legend:

``` r
pB <- ggplot(GCsbyCOG_Genome, aes(x = accessory_vs_core, y = num_corrected_genes, fill = fct_rev(Assignment))) +
  stat_summary(fun=sum ,geom="bar", position = "stack") +
  scale_x_discrete(labels = c("Soft/Core", "Accessory")) +
  scale_fill_manual(values = palette3) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 1500, by = 250)) +
  labs(fill=" ", x=" ", y= "Number of Gene Clusters") +
  theme_classic() +
  theme(axis.title = element_text(size = 9), axis.text = element_text(size=7), plot.margin=unit(c(10,0,10,20),"pt"), legend.key.size = unit(0.7, "line"), legend.text = element_text(size = 7), legend.box.margin = margin(10,20,10,10)) +
  guides(fill=guide_legend(ncol=1, title.position = "top", title.hjust = 0.5))
```

### Plots by Genome

Panel A in supplemental figure:

``` r
pC <- ggplot(filter(GCsbyCOG_Genome, accessory_vs_core == "Accessory"), aes(x=genome_name, y=num_corrected_genes, fill = fct_rev(COGs))) +
  stat_summary(fun=sum ,geom="bar", position = "stack") +
  scale_fill_manual(values = palette1) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(fill="COG Assignment", x="", y= "Number of Gene Clusters") +
  theme_classic() + 
  theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=8, angle=75, hjust=1)) +
  theme(legend.position = "none", plot.margin=unit(c(15,15,-10,20),"pt")) 
```

Panel B in supplemental figure:

``` r
pD <- ggplot(filter(GCsbyCOG_Genome %>% filter(COGs != "Uninformative", COGs !="Ambiguous", COGs != "Unclassified"), accessory_vs_core == "Accessory"), aes(x=genome_name, y=num_corrected_genes, fill = fct_rev(COGs))) +
  stat_summary(fun=sum ,geom="bar", position = "stack") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = palette2) + 
  labs(fill="COG Categories", x="", y= "Number of Informative Gene Clusters") +
  theme_classic() + 
  theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=8, angle=75, hjust=1)) +
  theme(legend.position="bottom", legend.key.size = unit(0.7, "line"), legend.text = element_text(size = 8), plot.margin=unit(c(0,15,0,20),"pt")) +
  guides(fill=guide_legend(ncol=2, title.position = "top", title.hjust = 0.5)) 
```

This plot is used for the grayscale legend:

``` r
pE <- ggplot(filter(GCsbyCOG_Genome, accessory_vs_core == "Accessory"), aes(x=genome_name, y=num_corrected_genes, fill = fct_rev(Assignment))) +
  stat_summary(fun=sum ,geom="bar", position = "stack") +
  scale_fill_manual(values = palette3) +
  scale_y_continuous(expand = c(0,0)) + 
  labs(fill="Accessory Genome COG Assignment", x="", y= "Number of Gene Clusters") +
  theme_classic() + 
  theme(axis.text.y = element_text(size=7), axis.text.x = element_text(size=8, angle=75, hjust=1)) +
  theme(legend.position="bottom", legend.key.size = unit(0.7, "line"), legend.text = element_text(size = 8), legend.title = element_text(face="bold", size = 12), plot.margin=unit(c(15,15,0,20),"pt")) +
  guides(fill=guide_legend(nrow=1, title.position = "top", title.hjust = -5))
```

This is used for the clade labels in supplemental figure:

``` r
pclades <- ggplot() +
  scale_y_continuous(limits = c(-1.5, 0.5), breaks = c(-1, 0)) +
  geom_segment(aes(x=0,xend=2.9,y=0,yend=0), color="#2c9b51ff") +
  geom_segment(aes(x=3,xend=4.9,y=0,yend=0), color="#a851a8ff") +
  geom_segment(aes(x=5,xend=9.9,y=0,yend=0), color="#e17139ff") +
  geom_segment(aes(x=10,xend=28,y=0,yend=0), color="#1b1d86ff") +
  annotate("text", x = 1.5, y = -1, label = "C1", fontface="bold", color="#2c9b51ff")+
  annotate("text", x = 4, y = -1, label = "C2", fontface="bold", color="#a851a8ff")+
  annotate("text", x = 7.3, y = -1, label = "C3", fontface="bold", color="#e17139ff")+
  annotate("text", x = 19, y = -1, label = "C4", fontface="bold", color="#1b1d86ff")+
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.margin=unit(c(0,0,5,20),"pt")) 
```

### Plots by COG Category

In order to represent the Soft/Core on the left of the plot with
absolute values per COG category we create `core.neg`; a negative
version of the `core` variable in GCsbyCOG\_CorevsAcc. Table converted
to the long format for plotting.

``` r
GCsbyCOG_CorevsAcc$core.neg <- -GCsbyCOG_CorevsAcc$Core
GCsbyCOG_CorevsAccLong <- gather(GCsbyCOG_CorevsAcc, accessory_vs_core, plotting, core.neg, Accessory)
```

Panel B in main figure:

``` r
pF <- ggplot(filter(GCsbyCOG_CorevsAccLong, COGs != "Uninformative", COGs != "Ambiguous", COGs != "Unclassified"), aes(x = COGs, y = plotting, fill = COGs)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = rev(palette2)) + 
  scale_x_discrete(position = "top") +
  labs(title= "COG Categories", x="", y= "Number of Gene Clusters") +
  coord_flip() +
  scale_y_continuous(limits = c(-170, 170), breaks = c(-150, -100, -50, 0, 50, 100, 150), label = c(150, 100, 50, 0, 50, 100, 150)) +
  geom_segment(aes(x=0,xend=19.5,y=0,yend=0)) +
  geom_label(aes(x = 20.5, y = -95, label = "      Soft/Core       "), fontface="bold", size=3, fill = "grey90", label.size=NA, label.padding = unit(0.3, "lines")) +
  geom_label(aes(x = 20.5, y = 95, label = "     Accessory      "), fontface="bold", size=3, fill = "grey90", label.size=NA, label.padding = unit(0.3, "lines")) +
  theme_classic() +
  theme(axis.title = element_text(size = 9), axis.text.x = element_text(size=7), axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "none", plot.margin=unit(c(5,10,10,25),"pt"), plot.title=element_text(face="bold", hjust=3, vjust=-3.9)) 

gpF <- ggplotGrob(pF)
```

    ## Warning: Removed 1 rows containing missing values (position_stack).

``` r
gpF$layout$clip[gpF$layout$name=="panel"] <- "off"
```

## Final Figures

Main figure:

``` r
pMain <- ggarrange(ggarrange(get_legend(pB), pA, ncol = 1, heights = c(0.2, 1)),
                   gpF, ncol = 2, labels = c("A", "B"), hjust=-0.5, vjust=2, widths = c(0.7, 2))

ggsave("Fig4_COG_summary.tiff", pMain, width = 9, height = 4, dpi = 150)
```

Supplemental Figure:

``` r
pSupple <- ggarrange(get_legend(pE),
                      ggarrange(pC+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), pD+theme(legend.position="none"), ncol = 1,  align = "v", labels = c("i", "ii"), hjust=-0.5, vjust=1, heights = c(1, 1)),
                      pclades, 
                      get_legend(pD), ncol = 1, heights = c(0.2, 2, 0.2, 0.6))

ggsave("FigS1D_COG_byGenome.tiff", pSupple, width = 8, height = 10, dpi = 150)
```
