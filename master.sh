#!/bin/bash

echo "Ejecutando proyecto..."

snakemake --use-conda --cores 10 --keep-going

echo "Creando la matriz de conteos..."

featureCounts -T 8 -a genomes/gut_db_clean.gff \
    -F GTF -t CDS -g locus_tag \
    -p -B -C -o align/counts/counts_matrix.tsv \
    align/bam/*.sorted.bam

echo "Ejecutando análisis de expresión diferencial y enriquecimiento funcional..."

Rscript rnaseq_enrichment_analysis

echo "¡Proyecto comletado!"