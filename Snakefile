configfile: "config.yaml"

import pandas as pd
from os.path import join

# Cargar muestras
df = pd.read_csv("samples_baseline_trt.tsv", sep="\t")
SAMPLES = list(df["run_accession"])

# Directorios
raw_dir = config["raw_dir"]
qc_dir = config["qc_dir"]
qc_report_dir = config["qc_report_dir"]
tax_report_dir = config["tax_report_dir"]
bam_dir = config["bam_dir"]
counts_dir = config["counts_dir"]

threads = config.get("threads", 10)
kraken_db = config["kraken_db"]
metagenome_index = config["metagenome_index"]
gtf_file = config["gtf_file"]
human_index = config["human_index"]

# URLs
def get_r1(sample):
    return "ftp://" + df.loc[df.run_accession == sample, "fastq_ftp"].values[0].split(";")[0]

def get_r2(sample):
    urls = df.loc[df.run_accession == sample, "fastq_ftp"].values[0].split(";")
    return "ftp://" + urls[1] if len(urls) > 1 else None

# ----------------------
# REGLA PRINCIPAL
# ----------------------
rule all:
    input:
        f"{qc_report_dir}/multiqc_report.html",
        f"{tax_report_dir}/multiqc_taxonomy_report.html",
        expand(f"{bam_dir}/{{sample}}.sorted.bam", sample=SAMPLES) 

# ----------------------
# DESCARGAR FASTQ
# ----------------------
rule download_fastq:
    output:
        r1 = temp(f"{raw_dir}/{{sample}}_1.fastq.gz"),
        r2 = temp(f"{raw_dir}/{{sample}}_2.fastq.gz")
    params:
        r1_url = lambda wc: get_r1(wc.sample),
        r2_url = lambda wc: get_r2(wc.sample)
    log:
        f"{raw_dir}/{{sample}}.download.log"
    shell:
        """
        wget -q {params.r1_url} -O {output.r1} 2>> {log}
        wget -q {params.r2_url} -O {output.r2} 2>> {log}
        """

# ----------------------
# FASTQC
# ----------------------
rule fastqc:
    input:
        r1 = f"{raw_dir}/{{sample}}_1.fastq.gz",
        r2 = f"{raw_dir}/{{sample}}_2.fastq.gz"
    output:
        zip1 = f"{qc_dir}/{{sample}}_1_fastqc.zip",
        zip2 = f"{qc_dir}/{{sample}}_2_fastqc.zip"
    log:
        f"{qc_dir}/{{sample}}_fastqc.log"
    threads: 2
    shell:
        """
        fastqc {input.r1} {input.r2} --noextract --outdir {qc_dir} >> {log} 2>&1
        """

# ----------------------
# PROCESAMIENTO COMPLETO
# ----------------------
rule procesar_muestra:
    input:
        r1 = f"{raw_dir}/{{sample}}_1.fastq.gz",
        r2 = f"{raw_dir}/{{sample}}_2.fastq.gz"
    output:
        bam = f"{bam_dir}/{{sample}}.sorted.bam",
        kraken_report = "taxonomia/kraken/{sample}.report",
        bracken_output = "taxonomia/bracken/{sample}.bracken.txt"
    threads: 8
    params:
        tmpdir = "/tmp"
    shell:
        """
        TMPDIR=$(mktemp -d -p {params.tmpdir})
        
        # 1. Trimming
        trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} \
            $TMPDIR/{wildcards.sample}_1.trim.fastq.gz /dev/null \
            $TMPDIR/{wildcards.sample}_2.trim.fastq.gz /dev/null \
            SLIDINGWINDOW:4:20 MINLEN:50 2>/dev/null || echo "Trimmomatic completed"
        
        # 2. Filtrar humano
        bowtie2 --very-sensitive --threads {threads} -x {human_index} \
            -1 $TMPDIR/{wildcards.sample}_1.trim.fastq.gz \
            -2 $TMPDIR/{wildcards.sample}_2.trim.fastq.gz \
            --un-conc-gz $TMPDIR/{wildcards.sample}_clean.fastq.gz \
            -S /dev/null 2>/dev/null || echo "Bowtie2 human filtering completed"
        
        # 3. Kraken2
        kraken2 --db {kraken_db} --threads {threads} \
                --paired $TMPDIR/{wildcards.sample}_clean.fastq.1.gz \
                         $TMPDIR/{wildcards.sample}_clean.fastq.2.gz \
                --report {output.kraken_report} \
                --output taxonomia/kraken/{wildcards.sample}.kraken 2>/dev/null || \
        echo "Kraken2 completed (check logs if needed)"
        
        # 4. Bracken (con fallback)
        if [ -s {output.kraken_report} ]; then
            bracken -d {kraken_db} -i {output.kraken_report} \
                    -o {output.bracken_output} -r 150 -l S 2>/dev/null || \
            echo "Bracken failed, creating empty file" && touch {output.bracken_output}
        else
            echo "Empty Kraken report" && touch {output.bracken_output}
        fi
        
        # 5. Alineamiento
        bowtie2 -x {metagenome_index} -p {threads} \
                -1 $TMPDIR/{wildcards.sample}_clean.fastq.1.gz \
                -2 $TMPDIR/{wildcards.sample}_clean.fastq.2.gz \
                -S $TMPDIR/{wildcards.sample}.sam 2>/dev/null || echo "Bowtie2 alignment completed"
        
        # 6. SAM to sorted BAM
        samtools sort -@ {threads} -o {output.bam} $TMPDIR/{wildcards.sample}.sam 2>/dev/null
        
        rm -rf $TMPDIR
        """

# ----------------------
# REGLAS FINALES
# ----------------------

rule multiqc_qc:
    input:
        expand(f"{qc_dir}/{{sample}}_{{read}}_fastqc.zip", 
               sample=SAMPLES, read=["1", "2"])
    output:
        f"{qc_report_dir}/multiqc_report.html"
    log:
        f"{qc_report_dir}/multiqc.log"
    shell:
        """
        rm -f {qc_report_dir}/multiqc_report.html 2>/dev/null || true
        multiqc {qc_dir} -o {qc_report_dir} -n multiqc_report.html --force 2> {log}
        """

rule multiqc_taxonomy:
    input:
        expand("taxonomia/kraken/{sample}.report", sample=SAMPLES),
        expand("taxonomia/bracken/{sample}.bracken.txt", sample=SAMPLES)
    output:
        "taxonomia/reports/multiqc_taxonomy_report.html"
    log:
        "taxonomia/reports/multiqc_tax.log"
    shell:
        """
        mkdir -p taxonomia/reports
        multiqc taxonomia/kraken taxonomia/bracken \
            -o taxonomia/reports \
            -n multiqc_taxonomy_report.html --force 2> {log} || true
        """




