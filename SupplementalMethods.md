*Dolosigranulum pigrum* Supplemental Methods
================

# SETUP

-   [Anvi’o v7,
    “hope”](https://github.com/merenlab/anvio/releases/tag/v7) was
    [installed](https://merenlab.org/2016/06/26/installation-v2/) in a
    Python environment called `anvio-7`.

-   [PPanGGOLiN v1.1.141](https://github.com/labgem/PPanGGOLiN/releases)
    was
    [installed](https://github.com/labgem/PPanGGOLiN/wiki/Installation)
    in a Python environment called `PPanGGOLiN.`

# ANVI’O GENOME CONTIGS DATABASES

## Reformating genome fasta files for Anvi’o

The fasta files need to be reformatted again for Anvi’o v7, since we
want to avoid a Config Error indicating the sequence “contains
characters that are not any of A, C, T, G, N, a, c, t, g, n.” Flag
–seq-type NT added to fix this error.

``` bash
#conda activate anvio-7
mkdir -p "analysis_Anvio7/Reformat_Dpig"
path_f="GENOMES/renamed"
path_o="analysis_Anvio7/Reformat_Dpig"

for file in $path_f/*.f*; do
    FILENAME=`basename ${file%.*}`
    anvi-script-reformat-fasta -o $path_o/$FILENAME.fa --min-len 0 --simplify-names $file --seq-type NT; 
done
```

## Parsing .gff files

The next step is to parse Prokka annotated genomes to add external gene
calls and functions into Anvi’o. The input is the annotation in GFF3
format and outputs are tab-delimited text files, one for gene calls and
one for annotations. This is done with the script `gff_parser.py`.

The parsing of the Prokka genomes needs to be done again using the new
version of the `gff_parser.py` script in order to generate the right
columns expected in the genome calls files.

``` bash
#conda activate anvio-7
mkdir -p "analysis_Anvio7/Parsed_prokka_Dpig"
path_f="GENOMES/Prokka_Out_Renamed"
path_o="analysis_Anvio7/Parsed_prokka_Dpig"

for file in $path_f/*.gff; do
    FILENAME=`basename ${file%.*}`
    python analysis_Anvio7/gff_parser.py --gene-calls $path_o/calls_$FILENAME.txt --annotation $path_o/annot_$FILENAME.txt $file;
done
```

## Generating contigs databases

In this step the gene calls from Prokka and the reformatted fasta files
get imported to generate Anvi’o contigs databases. In the past (Coryne
genomes) we got a lot of early stop codon errors. Therefore, we add the
–ignore-internal-stop-codons flag. Based on the terminal output
(`anvi-gen-contigs-database_log.txt`), all genomes have 1-2 genes with
internal stop codons:

``` bash
#conda activate anvio-7
mkdir -p "analysis_Anvio7/Contigs_db_prokka_Dpig"
path_f="analysis_Anvio7/Reformat_Dpig"
path_o="analysis_Anvio7/Contigs_db_prokka_Dpig"
path_e="analysis_Anvio7/Parsed_prokka_Dpig"

for file in $path_f/*.fa; do
    FILENAME=`basename ${file%.*}`
    anvi-gen-contigs-database -f $file \
                              -o $path_o/$FILENAME.db \
                              --external-gene-calls $path_e/calls_prokka_$FILENAME.txt \
                              --ignore-internal-stop-codons \
                              -n $FILENAME;
done
```

## Importing functional annotation

### Prokka annotation

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"
path_e="analysis_Anvio7/Parsed_prokka_Dpig"

for file in $path_f/*.db; do
    FILENAME=`basename ${file%.*}`
    anvi-import-functions -c $file \
                          -i $path_e/annot_prokka_$FILENAME.txt
      
done
```

### COG annotation

First you need to download and set up the NCBI’s Clusters of Orthologous
Groups database using `anvi-setup-ncbi-cogs`: this is something you will
do only once.

Now the `anvi-run-ncbi-cogs` command can be used to annotate the .db
genome files against the NCBI’s COGs database: It uses DIAMOND “fast” as
default, but it is recommended to use the “sensitive” option. This will
take considerably much longer. Blastp can also be used instead of
DIAMOND.

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-ncbi-cogs -T 4 --sensitive -c $file;
done
```

### KEGG annotation

First you need to download and set up the KEGG KOfam database using
`anvi-setup-kegg-kofams`: this is something you will do only once.

The program `anvi-run-kegg-kofams` annotates genes in a given anvi’o
contigs database with KEGG Orthology (KO) numbers (via hits to the KEGG
KOfam database).

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-kegg-kofams -T 4 -c $file;
done
```

### PFAM annotation

First you need to download and set up the PFAM database using
`anvi-setup-pfams`: this is something you will do only once.

The program `anvi-run-kegg-pfam` annotates genes in a given anvi’o
contigs database with the PFAM database.

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-pfams -T 4 -c $file;
done
```

### Running HMMs

When you run the following command the occurrence of bacterial
single-copy genes across your contigs is all added to your contigs
database.

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-hmms -T 4 -c $file;
done
```

### Annotation Stats

``` bash
#conda activate anvio-7
mkdir -p "analysis_Anvio7/Contigs_stats"

path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"
path_o="analysis_Anvio7/Contigs_stats"

for file in $path_f/*.db; do
    FILENAME=`basename ${file%.*}`
    anvi-display-contigs-stats $file --report-as-text -o $path_o/stats_$FILENAME.txt;
done
```

# ANVI’O PANGENOME ANALYSIS

## Generating an Anvi’o genomes storage

The program `anvi-gen-genomes-storage` requires you to provide names and
locations of the genomes to be included in this storage in a .txt file.
We create the AnvioGenomesSummary003.txt file with the renamed file
names and paths.

**IMPORTANT:** The next command generates big files as output that are
too big to be sync using GitHub; please don’t commit to GitHub.

``` bash
#conda activate anvio-7
anvi-gen-genomes-storage -e analysis_Anvio7/Contigs_db_prokka_Dpig/AnvioGenomesSummary003.txt -o analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db --gene-caller Prodigal
```

## Running the pangenomic analysis

Once the genomes storage is ready, you can use the program
`anvi-pan-genome` to run the actual pangenomic analysis. There are many
parameters that can be adjusted for this step. The most relevant are:

-   Basic flags, file management:

    -   `-g` your genome storage location
    -   `-o` your output files location
    -   `-n` is the name for your project
    -   `--num-threads` use double the number of cores on your computer
        (for a dual-cored laptop use 4)

-   Advanced flags, algorithm parameters:

    -   `--use-ncbi-blast` Meren strongly suggest to use blastp instead
        of DIAMOND (the default is DIAMONDfast)
    -   `--mcl-inflation` 0 to 10, the default is 2. Meren recommends
        using 2 for distantly related genomes (i.e., genomes classify
        into different families or farther), and 10 for very closely
        related genomes (i.e., ‘strains’ of the same ‘species’)
    -   `--minbit` 0 to 1, the default is 0.5. The minbit heuristic
        provides a mean to eliminate weak matches between two amino acid
        sequences. It seems to take into account both similarity and
        coverage
    -   `--min-occurrence` the default is 1. Increase if you want to
        ignore singletons, doubletons… or define a “Core” pangenome

**IMPORTANT:** The next command generates several intermediate files
that are too big to be sync using GitHub; please don’t commit those
files to GitHub.

``` bash
#conda activate anvio-7
anvi-pan-genome -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db -o analysis_Anvio7/Pangenomic_Results_Dpig/ -n PAN_DPIG_prokka --use-ncbi-blast --mcl-inflation 10  --num-threads 4
```

## Displaying the pangenome

Start an anvi’o server to display a pan-genome. This will display the
pangenome in your default internet browser.

`-p` input for PAN.db file `-g` genomes storage file relating to that
PAN.db

``` bash
#conda activate anvio-7
sudo anvi-display-pan -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db
```

Once it opens in the browser, you can adjust image settings and click
the “DRAW” icon to view it. For help using the interactive interface see
<http://merenlab.org/2016/02/27/the-anvio-interactive-interface/>

## Creating bin collections (Core vs Accessory)

In order to define Core vs Accessory pangenome we use the interactive
interface. In the “Bins” tab we can create and name as many bins as we
want and store then in a bin collection. In the “Search” tab I used
“Search gene clusters using filters” and “Append splits to selected bin”
to create the following bins:

-   Core: 1298/2905 (44.7%) gene clusters in all 28 genomes

    -   Min number of genomes gene cluster occurs = 28

-   SC Core: 1111/2905 (38.2%) gene clusters detected in single copy in
    all 28 genomes

    -   Min number of genomes gene cluster occurs = 28
    -   Max number of genes from each genome = 1

-   MC Core: 1298-1111=187 (6.4%) detected in multiple copy in all 28
    genomes

-   SoftCore: 90/2905 (3.1%) gene clusters detected in 26-27 genomes

    -   Min number of genomes gene cluster occurs = 26
    -   Max number of genomes gene cluster occurs = 27

-   Shell: 805/2905 (27.7%) gene clusters detected in 3-25 genomes

    -   Min number of genomes gene cluster occurs = 3
    -   Max number of genomes gene cluster occurs = 25

-   Cloud: 712/2905 (24.5%) gene clusters detected in only 2 genomes or
    less

    -   Max number of genomes gene cluster occurs = 2

We stored the 5 bins in a new collection named **CorevsAccessory**

We can see the collections and bins created in a pangenomic database
using the command `anvi-show-collections-and-bins`

``` bash
#conda activate anvio-7
anvi-show-collections-and-bins --debug -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db
```

## Summarizing the pangenome

We can summarized the pangenome using the program `anvi-summarize`:

``` bash
#conda activate anvio-7
anvi-summarize -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db \
               -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db \
               -C CorevsAccessory \
               -o analysis_Anvio7/Pangenomic_Results_Dpig/Dpig-PAN-SUMMARY
```

The resulting summary folder contains the file
PAN\_DPIG\_prokka\_gene\_clusters\_summary.txt.gz that links each gene
to gene clusters, genomes, functions, and bins selected from the
interface.

## Displaying functions

``` bash
anvi-db-info analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db
```

gene\_function\_sources ……………………:
KOfam,Pfam,KEGG\_Module,KEGG\_Class,COG20\_PATHWAY,COG20\_FUNCTION,COG20\_CATEGORY,Prokka:Prodigal

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-db-info $file;
done
```

# PPanGGOLiN-Anvi’o 7 PIPELINE

## Importing Prokka annotations and Anvi’o 7 clusters

In order to import the Anvi’o clustering into PPanGGOLiN we need:

1.  A .tsv file listing in the first column the gene family names, and
    in the second column the gene ID that is used in the annotation
    files.
2.  The annotated genomes with gene IDs that match the ones listed in
    the previous .tsv file.

Using `anvi-summarize` we can generate files that link each gene to gene
clusters, genomes, functions, and bins selected from the interface. And
we can subset the info needed to generate the .tsv file from it:

``` r
Dpig_Anvio7 <- read_delim("analysis_Anvio7/Pangenomic_Results_Dpig/Dpig-PAN-SUMMARY/PAN_DPIG_prokka_gene_clusters_summary.txt.gz", "\t")
Dpig_Anvio7 <- Dpig_Anvio7 %>% 
  unite(new_id, genome_name:gene_callers_id, sep = "___", remove = FALSE)
Clusters_Dpig_Anvio7 <- select(Dpig_Anvio7, c(gene_cluster_id, new_id))
```

``` r
write_delim(Clusters_Dpig_Anvio7, "analysis_PPanGGOLiN_Anvio7/Clusters_Dpig_Anvio7.tsv", col_names=FALSE)
```

**IMPORTANT:** The gene\_callers\_id provided by Anvi’o don’t match the
original ones in the Prokka annotation.

The Prokka annotated genomes were parsed into two text files, one for
gene calls and one for annotations, with the script `gff_parser.py`. By
default Prokka annotates also tRNAs, rRNAs and CRISPR regions. However,
`gff_parser.py` will only utilize open reading frames reported by
Prodigal in the Prokka output in order to be compatible with the
pangenomic Anvi’o pipeline. While parsing new gene\_callers\_id are
generated only for the ORFs that will be imported into Anvi’o.
Fortunately `anvi-get-sequences-for-gene-calls` can be used to export
new .gff files with only the ORFs and these can be used for PPanGGOLiN,
instead of the original Prokka ones:

``` bash
#conda activate anvio-7
mkdir -p "analysis_Anvio7/Exported_gffs"

path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"
path_o="analysis_Anvio7/Exported_gffs"

for file in $path_f/*.db; do
    FILENAME=`basename ${file%.*}`
    anvi-get-sequences-for-gene-calls -c $file --export-gff3 \
                                      -o $path_o/Anvio7_$FILENAME.gff
      
done
```

I create lists with the files names for both the Anvi’o exported .gff
and the renamed .fasta files and use them to run the `annotate`
subcommand. For the `cluster` subcommand I use the .tsv file created
before from the `anvi-summarize` output.

``` bash
#conda activate PPanGGOLiN
ppanggolin annotate --anno analysis_PPanGGOLiN_Anvio7/Anvio7GenomesExported.gff.txt --fasta analysis_PPanGGOLiN_Anvio7/Anvio7GenomesReformatted.fa.txt -o analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7 --basename FromAnvio7

ppanggolin cluster -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5 --clusters analysis_PPanGGOLiN_Anvio7/Clusters_Dpig_Anvio7.tsv --infer_singletons
```

## Graphing and Partitioning

The `graph` subcommand has only a single other option, which is ‘-r’ or
‘–remove\_high\_copy\_number’. If used, it will remove the gene families
that are too duplicated in your genomes. This is useful if you want to
visualize your pangenome afterward and want to remove the biggest hubs
to have a clearer view. It can also be used to limit the influence of
very duplicated genes such as transposase or ABC transporters in the
partition step.

We let the `partition` subcommand statistical criterion find the optimal
number of partitions.

``` bash
#conda activate PPanGGOLiN
ppanggolin graph -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5

ppanggolin partition -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5
```

## Writing outpus

For details on PPanGGOLiN outputs see:
<https://github.com/labgem/PPanGGOLiN/wiki/Outputs>

``` bash
#conda activate PPanGGOLiN
ppanggolin write -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5 -o analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7 --light_gexf --gexf  --csv --Rtab --stats --partitions --projection --families_tsv -f
```

## Finding Regions of Genome Plasticity

For details on RGPs see:
<https://github.com/labgem/PPanGGOLiN/wiki/Regions-of-Genome-Plasticity>

``` bash
#conda activate PPanGGOLiN
ppanggolin rgp -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5
ppanggolin spot -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5 --draw_hotspots -o analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/spots

ppanggolin write -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5 -o analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7 --regions --spots -f
```

## Summary of used parameters

The `info` subcommand indicates, for each steps of the analysis, the
PPanGGOLiN parameters that were used and the source of the data if
appropriate.

``` bash
#conda activate PPanGGOLiN
ppanggolin info -p analysis_PPanGGOLiN_Anvio7/OutputFromAnvio7/FromAnvio7.h5 --parameters
```

The output for this command was:

-   annotation

    -   read\_annotations\_from\_file : True

-   cluster

    -   read\_clustering\_from\_file : True
    -   infer\_singletons : True

-   graph

    -   removed\_high\_copy\_number\_families : False

-   partition

    -   beta : 2.5
    -   free\_dispersion : False
    -   max\_node\_degree\_for\_smoothing : 10
    -   computed\_K : True
    -   K : 3

-   RGP

    -   persistent\_penalty : 3
    -   variable\_gain : 1
    -   min\_length : 3000
    -   min\_score : 4
    -   dup\_margin : 0.05

-   spots

    -   set\_size : 3
    -   overlapping\_match : 2
    -   exact\_match : 1
