# RNA-seq Analysis Pipeline: B. longum AH1206
## Pipeline para procesar datos de metatranscriptómica (PRJNA324129) usando Snakemake y R.

### Requisitos:
- Conda / Mamba

### Bases de Datos (Configurar en config.yaml):

- MiniKraken2 DB or other
- GRCh38 Indexed Human Genome

### Instalación:

```{bash}
conda env create -f envs/bioinformatics.yml
conda activate bioinfo_probiotic
```

### Configuración Inicial
Ejecutar el script de setup para descargar genomas de referencia, crear carpetas, indexar con Bowtie2 y preparar los scripts en la raíz:

```{bash}
chmod +x setup.sh
bash setup.sh
```

### Ejecución Completa
El script `master.sh` lanza todo el proceso de forma secuencial:

1) Snakemake: QC, trimming, filtrado de host , análisis taxonómico y mapeo.

2) featureCounts: Generación de la matriz de conteos.

3) Rscript: Análisis de expresión diferencial y enriquecimiento.

```{bash}
chmod +x master.sh
bash master.sh
```

Tras la ejecución del pipeline, podrá consultar por pantalla las tablas HTML y gráficos de R, y observar los análisis MultiQC de calidad (disponible en `control/calidad/reports/`) y de taxonomía (disponible en `taxonomia/reports/`). 

### Notas

Este pipeline analiza 88 muetras, consultables en `samples_baseline_trt.tsv` y puede resultar lento de ejecutar. Considere en usar la opción `snakemake --cores 10 --rerun-incomplete` de Snakemake si no obtiene un resultado satisfactorio del pipeline.



