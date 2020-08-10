## Generating Anvi’o genomes storage
conda activate anvio-6.2

cd /YOUR/PATH/TO/AnvioPangenome 

anvi-gen-genomes-storage -e GenomeDatabases/Dpi_external_genomes.txt \
                         -o PROKKA_GENOMES.db --gene-caller Prodigal


##Use this to display the pangenome on Anvi'o##

anvi-display-pan -g PROKKA_GENOMES.db \
                 -p Dpi_Prokka_Pan_t28-PAN.db
                 