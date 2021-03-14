"Genomic Stability and Genetic Defense Systems in *Dolosigranulum
pigrum* a Candidate Beneficial Bacterium from the Human Microbiome"
Supplemental Methods
================

# ANVI’O PANGENOME ANALYSIS

[Anvi’o v7, “hope”](https://github.com/merenlab/anvio/releases/tag/v7)
([Eren et al. 2015](#ref-eren2015), [2020](#ref-eren2020)) was
[installed](https://merenlab.org/2016/06/26/installation-v2/) in a
Python environment called `anvio-7`.

We followed a pipeline modified from the tutorial: [Anvi’o workflow for
microbial
pangenomics](https://merenlab.org/2016/11/08/pangenomics-v2/ "An anvi'o workflow for microbial pangenomics").

## Reformating genome fasta files

The NCBI submitted .gb genome files were converted to .fa using SnapGene
and renamed as indicated on `GENOMES/GenomeNames.txt`. The renamed
genome files are available in the `GENOMES/renamed` folder.

These .fa genomes files were reformatted to be compatible with Anvi’o
v7. Flag `--seq-type NT` was added to avoid a Config Error indicating
the sequence “contains characters that are not any of A, C, T, G, N, a,
c, t, g, n.”

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

## Importing Prokka annotations into Anvi’o

The Prokka annotation was performed with [Prokka
v1.14.6](https://github.com/tseemann/prokka/releases/tag/v1.14.6)
([Seemann 2014](#ref-seemann2014)) with default parameters, including
gene recognition and translation initiation site identification with
Prodigal ([Hyatt et al. 2010](#ref-hyatt2010)). Prokka annotations were
imported into Anvi’o following the recommended [Anvi’o pipeline for
Prokka annotated
genomes](https://merenlab.org/2017/05/18/working-with-prokka/) described
in the next steps.

### Parsing .gff files

The Prokka annotated .gff files are parsed into two tab-delimited text
files; one containing the Prodigal gene calls data and another with the
functional annotations. This is done with the script `gff_parser.py`.

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

### Generating contigs databases

In this step the gene calls from Prokka and the reformatted fasta files
get imported to generate Anvi’o contigs databases. We add the
`--ignore-internal-stop-codons` flag to avoid errors downstream the
pipeline. Only 1-2 genes with internal stop codons were identified in
each genome.

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

### Importing Prokka annotation

The functional annotation part of the parsed Prokka files gets added to
the Anvi’o contigs databases (.db files)

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

## Importing other functional annotations

### COG annotation

The `anvi-run-ncbi-cogs` command was used to annotate the .db genome
files against the NCBI’s COGs database 2020 release ([Tatusov
1997](#ref-tatusov1997); [Galperin et al. 2020](#ref-galperin2020)). We
used the DIAMOND ([Buchfink, Xie, and Huson 2014](#ref-buchfink2014))
“sensitive” option instead of the default “fast” one.

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-ncbi-cogs -T 4 --sensitive -c $file;
done
```

### KEGG annotation

The program `anvi-run-kegg-kofams` annotates the .db genomes with KEGG
Orthology (KO) numbers via hits to the KEGG KOfam database ([M. Kanehisa
2000](#ref-Kanehisa2000); [Minoru Kanehisa et al.
2015](#ref-Kanehisa2016)).

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-kegg-kofams -T 4 -c $file;
done
```

### PFAM annotation

The program `anvi-run-pfams` annotates the .db genomes with the PFAM
database ([Mistry et al. 2020](#ref-mistry2020)).

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-pfams -T 4 -c $file;
done
```

### Running HMMs

The program
[`anvi-run-hmms`](https://merenlab.org/software/anvio/help/main/programs/anvi-run-hmms/)
stores hmm-hits into the .db genomes files ([Eddy 2011](#ref-eddy2011)).

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-run-hmms -T 4 -c $file;
done
```

### Displaying functions

The program
[`anvi-db-info`](https://merenlab.org/software/anvio/help/main/programs/anvi-db-info/)
displays information about an Anvi’o database. The output of this
command included in `analysis_Anvio7/AnvioGenomesSummaryFuntions.txt`
has the details of all the functional annotations included in the .db
files.

``` bash
#conda activate anvio-7
path_f="analysis_Anvio7/Contigs_db_prokka_Dpig"

for file in $path_f/*.db; do
    anvi-db-info $file;
done
```

## Generating an Anvi’o genomes storage

The program `anvi-gen-genomes-storage` requires you to provide names and
locations of the genomes to be included in a genome storage in a .txt
file. The `analysis_Anvio7/AnvioGenomesSummary003.txt` file includes the
file names and file paths.

``` bash
#conda activate anvio-7
anvi-gen-genomes-storage -e analysis_Anvio7/Contigs_db_prokka_Dpig/AnvioGenomesSummary003.txt -o analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db --gene-caller Prodigal
```

## Running the pangenomic analysis

The program `anvi-pan-genome`:

1.  Calculates similarities of each amino acid sequence in every genome
    against every other amino acid sequence. We used the more sensitive
    blastp search (`--use-ncbi-blast`) ([Altschul et al.
    1990](#ref-altschul1990)) instead of the default DIAMOND fast
    option.

2.  Uses ‘minbit heuristic’ ([Benedict et al. 2014](#ref-benedict2014))
    to filter weak hits based on the aligned fraction between the two
    reads.

3.  Uses the MCL algorithm ([van Dongen and Abreu-Goodger
    2011](#ref-vandongen2011)) to identify gene clusters in the
    remaining blastp search results. `--mcl-inflation` was set to 10, as
    recommended for very closely related genomes (i.e., ‘strains’ of the
    same ‘species’).

``` bash
#conda activate anvio-7
anvi-pan-genome -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db -o analysis_Anvio7/Pangenomic_Results_Dpig/ -n PAN_DPIG_prokka --use-ncbi-blast --mcl-inflation 10  --num-threads 4
```

## Displaying the pangenome

This will display the pangenome in your default internet browser. Once
it opens in the browser, you can adjust image settings and click the
“DRAW” icon to view it. For help using the interactive interface see
<http://merenlab.org/2016/02/27/the-anvio-interactive-interface/>

``` bash
#conda activate anvio-7
sudo anvi-display-pan -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db
```

## Creating bin collections (Core vs Accessory)

In order to define Core vs Accessory pangenome we use Anvi’o interactive
interface. In the “Bins” tab we can create and name as many bins as we
want and store then in a bin collection. In the “Search” tab we used
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

We summarized the pangenome using the program `anvi-summarize`:

``` bash
#conda activate anvio-7
anvi-summarize -p analysis_Anvio7/Pangenomic_Results_Dpig/PAN_DPIG_prokka-PAN.db \
               -g analysis_Anvio7/Pangenomic_Results_Dpig/Anvio005-GENOMES.db \
               -C CorevsAccessory \
               -o analysis_Anvio7/Pangenomic_Results_Dpig/Dpig-PAN-SUMMARY
```

The resulting summary folder contains the file
`analysis_Anvio7/Pangenomic_Results_Dpig/Dpig-PAN-SUMMARY/PAN_DPIG_prokka_gene_clusters_summary.txt.gz`
that links each gene to gene clusters, genomes, functions, and bins
selected from the interface.

# <u>REFERENCES</u>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-altschul1990" class="csl-entry">

Altschul, Stephen F., Warren Gish, Webb Miller, Eugene W. Myers, and
David J. Lipman. 1990. “Basic Local Alignment Search Tool.” *Journal of
Molecular Biology* 215 (3): 403–10.
<https://doi.org/10.1016/s0022-2836(05)80360-2>.

</div>

<div id="ref-benedict2014" class="csl-entry">

Benedict, Matthew N, James R Henriksen, William W Metcalf, Rachel J
Whitaker, and Nathan D Price. 2014. “ITEP: An Integrated Toolkit for
Exploration of Microbial Pan-Genomes.” *BMC Genomics* 15 (1): 8.
<https://doi.org/10.1186/1471-2164-15-8>.

</div>

<div id="ref-buchfink2014" class="csl-entry">

Buchfink, Benjamin, Chao Xie, and Daniel H Huson. 2014. “Fast and
Sensitive Protein Alignment Using DIAMOND.” *Nature Methods* 12 (1):
59–60. <https://doi.org/10.1038/nmeth.3176>.

</div>

<div id="ref-eddy2011" class="csl-entry">

Eddy, Sean R. 2011. “Accelerated Profile HMM Searches.” Edited by
William R. Pearson. *PLoS Computational Biology* 7 (10): e1002195.
<https://doi.org/10.1371/journal.pcbi.1002195>.

</div>

<div id="ref-eren2015" class="csl-entry">

Eren, A. Murat, Özcan C. Esen, Christopher Quince, Joseph H. Vineis,
Hilary G. Morrison, Mitchell L. Sogin, and Tom O. Delmont. 2015.
“Anvi’o: An Advanced Analysis and Visualization Platform for ‘Omics
Data.” *PeerJ* 3 (October): e1319. <https://doi.org/10.7717/peerj.1319>.

</div>

<div id="ref-eren2020" class="csl-entry">

Eren, A. Murat, Evan Kiefl, Alon Shaiber, Iva Veseli, Samuel E. Miller,
Matthew S. Schechter, Isaac Fink, et al. 2020. “Community-Led,
Integrated, Reproducible Multi-Omics with Anvi’o.” *Nature Microbiology*
6 (1): 3–6. <https://doi.org/10.1038/s41564-020-00834-3>.

</div>

<div id="ref-galperin2020" class="csl-entry">

Galperin, Michael Y, Yuri I Wolf, Kira S Makarova, Roberto Vera Alvarez,
David Landsman, and Eugene V Koonin. 2020. “COG Database Update: Focus
on Microbial Diversity, Model Organisms, and Widespread Pathogens.”
*Nucleic Acids Research* 49 (D1): D274–81.
<https://doi.org/10.1093/nar/gkaa1018>.

</div>

<div id="ref-hyatt2010" class="csl-entry">

Hyatt, Doug, Gwo-Liang Chen, Philip F LoCascio, Miriam L Land, Frank W
Larimer, and Loren J Hauser. 2010. “Prodigal: Prokaryotic Gene
Recognition and Translation Initiation Site Identification.” *BMC
Bioinformatics* 11 (1). <https://doi.org/10.1186/1471-2105-11-119>.

</div>

<div id="ref-Kanehisa2000" class="csl-entry">

Kanehisa, M. 2000. “KEGG: Kyoto Encyclopedia of Genes and Genomes.”
*Nucleic Acids Research* 28 (1): 27–30.
<https://doi.org/10.1093/nar/28.1.27>.

</div>

<div id="ref-Kanehisa2016" class="csl-entry">

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

<div id="ref-seemann2014" class="csl-entry">

Seemann, T. 2014. “Prokka: Rapid Prokaryotic Genome Annotation.”
*Bioinformatics* 30 (14): 2068–69.
<https://doi.org/10.1093/bioinformatics/btu153>.

</div>

<div id="ref-tatusov1997" class="csl-entry">

Tatusov, R. L. 1997. “A Genomic Perspective on Protein Families.”
*Science* 278 (5338): 631–37.
<https://doi.org/10.1126/science.278.5338.631>.

</div>

<div id="ref-vandongen2011" class="csl-entry">

van Dongen, Stijn, and Cei Abreu-Goodger. 2011. “Using MCL to Extract
Clusters from Networks.” In, 281–95. Springer New York.
<https://doi.org/10.1007/978-1-61779-361-5_15>.

</div>

</div>
