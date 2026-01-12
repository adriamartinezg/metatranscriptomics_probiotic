#!/bin/bash

echo "Preparando el escenario de carpetas..."

mkdir -p raw trim clean qc
mkdir -p control_calidad/reports
mkdir -p taxonomia/kraken taxonomia/bracken taxonomia/reports
mkdir -p align/bam align/counts genomes

echo "Estructura de directorios creada correctamente."


datasets download genome taxon Lactobacillus --reference --assembly-level complete --include gff3,genome --filename lacto.zip
datasets download genome taxon Bifidobacterium --reference --assembly-level complete --include gff3,genome --filename bifido.zip
datasets download genome taxon Bacteroides --reference --assembly-level complete --include gff3,genome --filename bacte.zip
datasets download genome taxon Faecalibacterium --reference --assembly-level complete --include gff3,genome --filename faecali.zip
datasets download genome taxon Escherichia --reference --assembly-level complete --include gff3,genome --filename escher.zip

unzip -o "*.zip"

find ncbi_dataset -name "*.fna" -exec cat {} + > genomes/gut_db.fna
find ncbi_dataset -name "*.gff" -exec grep -v "^#" {} + > genomes/gut_db.gff

bowtie2-build genomes/gut_db.fna genomes/gut_db_idx

sed 's/^.*genomic.gff://' genomes/gut_db.gff > genomes/gut_db_temp.gff
grep -v "^#" genomes/gut_db_temp.gff | grep -E "\t(CDS|gene|region)\t" > genomes/gut_db_clean.gff
rm genomes/gut_db_temp.gff

rm -rf ncbi_dataset *.zip

if [ -f "scripts/Snakefile" ]; then
    mv scripts/Snakefile .
    echo "Snakefile movido a la raíz."
fi

if [ -f "scripts/rnaseq_enrichment_analysis.R" ]; then
    mv scripts/rnaseq_enrichment.R .
    echo "Script de R movido a la raíz."
fi

# Instalación paquete Tami
Rscript -e "if (!require('tami')) { install.packages('http://www.uv.es/ayala/software/tami_1.0.tar.gz', repos=NULL, type='source') }"