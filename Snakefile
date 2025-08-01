SAMPLES = ["SRR23609077", "SRR23609078", "SRR23609079", "SRR23609080", "SRR23609081", "SRR23609082", "SRR23609083", "SRR23609084", "SRR23609085", "SRR23609086"]
READS = ["1", "2"]

def get_fastqc_files(wildcards):
    return expand("output/fastqc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READS)

def get_trimmed_fastqc_files(wildcards):
    return expand("output/fastqc/{sample}_{read}_trimmed_fastqc.zip", sample=SAMPLES, read=READS)


rule all:
    input:
        # FastQC na surowych odczytach
        expand("output/fastqc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READS),
        expand("output/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=READS),

        # MultiQC na surowych odczytach
        "output/multiqc/raw/raw.html",

        # Trimming
        expand("output/trimmed/{sample}_1_trimmed.fastq.gz", sample=SAMPLES, read=READS),
        expand("output/trimmed/{sample}_2_trimmed.fastq.gz", sample=SAMPLES, read=READS),

        # FastQC na przycietych odczytach
        expand("output/fastqc/{sample}_{read}_trimmed_fastqc.zip", sample=SAMPLES, read=READS),
        expand("output/fastqc/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=READS),

        # MultiQC na przycietych odczytach
        "output/multiqc/trimmed/trimmed.html",

        # Reference genome index
        "reference_genome/NC_045512_2.fasta.bwt",

        # Mapowanie
        expand("output/mapped/{sample}.sam", sample=SAMPLES),
        expand("output/mapped/{sample}.bam", sample=SAMPLES),
        expand("output/sorted/{sample}_sorted.bam", sample=SAMPLES),
        expand("output/sorted/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand("output/alignment_score/{sample}_sorted_alignment.txt", sample=SAMPLES),
        expand("output/variants/{sample}.vcf", sample=SAMPLES)

rule variant_calling:
    input:
        sorted_file="output/sorted/{sample}_sorted.bam",
        ref_file="reference_genome/NC_045512_2.fasta"
    output:
        variants_file="output/variants/{sample}.vcf"
    container:
            "docker://staphb/bcftools"
    shell:
        """
        mkdir -p output/variants
        bcftools mpileup -Ou -f {input.ref_file} {input.sorted_file} | \
        bcftools call -mv -Ov -o {output.variants_file} 
        """

rule alignment_score:
    input:
        sorted_file="output/sorted/{sample}_sorted.bam"
    output:
        score_file="output/alignment_score/{sample}_sorted_alignment.txt"
    container:
        "docker://staphb/samtools"
    shell:
        """
        mkdir -p output/alignment_score
        samtools coverage {input.sorted_file} > {output.score_file}
        """
       
rule sorting:
    input:
        mapped_file="output/mapped/{sample}.bam"
    output:
        sorted_file="output/sorted/{sample}_sorted.bam",
        sorted_indexed_file="output/sorted/{sample}_sorted.bam.bai"
    container:
        "docker://staphb/samtools"
    shell:
        """
        mkdir -p output/sorted/indexed
        samtools sort -o {output.sorted_file} {input.mapped_file} 
        samtools index {output.sorted_file}
        """ 

rule indexing:
    input:
        ref_file="reference_genome/NC_045512_2.fasta"
    output:
        index_ref_file="reference_genome/NC_045512_2.fasta.bwt"
    container:
        "docker://staphb/bwa"
    shell:
        """
        bwa index {input.ref_file}
        """

rule converting_sam:
    input:
        mapped_file="output/mapped/{sample}.sam"
    output:
        mapped_file_bam="output/mapped/{sample}.bam"
    container:
        "docker://staphb/samtools"
    shell:
        """
        samtools view -bS {input.mapped_file} > {output.mapped_file_bam} 
        """

rule mapping:
    input: 
        trimmed_file1="output/trimmed/{sample}_1_trimmed.fastq.gz",
        trimmed_file2="output/trimmed/{sample}_2_trimmed.fastq.gz",
        index_ref_file="reference_genome/NC_045512_2.fasta"
    output:
        mapped_file="output/mapped/{sample}.sam"
    container:
        "docker://staphb/bwa"
    shell:
        """
        mkdir -p output/mapped
        bwa mem {input.index_ref_file} {input.trimmed_file1} {input.trimmed_file2} > {output.mapped_file}
        """

rule trimmomatic:
    input:
        raw_file1="input/{sample}_1.fastq.gz",
        raw_file2="input/{sample}_2.fastq.gz"
    output:
        trimmed_file1="output/trimmed/{sample}_1_trimmed.fastq.gz",
        trimmed_file2="output/trimmed/{sample}_2_trimmed.fastq.gz"
    container:
        "docker://staphb/trimmomatic:0.38"
    shell:
        """
        mkdir -p output/trimmed
        trimmomatic PE -threads 4 -phred33 \
        {input.raw_file1} {input.raw_file2} \
        {output.trimmed_file1} /dev/null \
        {output.trimmed_file2} /dev/null \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule multiqc_raw:
    input:
        get_fastqc_files
    output:
        "output/multiqc/raw/raw.html"
    container:
        "docker://multiqc/multiqc:dev"
    shell:
        "multiqc output/fastqc --filename multiqc_report.html -o output/multiqc/raw -f -n raw"

rule multiqc_trimmed:
    input:
        get_trimmed_fastqc_files
    output:
        "output/multiqc/trimmed/multiqc_report.html"
    container:
        "docker://multiqc/multiqc:dev"
    shell:
        "multiqc output/fastqc --filename multiqc_report.html -o output/multiqc/trimmed -f -n trimmed"

rule fastqc:
    input:
        fastq_file="input/{sample}_{read}.fastq.gz"
    output:
        zip_files="output/fastqc/{sample}_{read}_fastqc.zip",
        html_files="output/fastqc/{sample}_{read}_fastqc.html"
    threads: lambda wildcards, attempt: min(attempt * 2, 8)
    log:
        "output/fastqc/{sample}_{read}.fastqc.log"
    container:
        "docker://staphb/fastqc"
    shell:
        "fastqc -t {threads} --outdir output/fastqc {input.fastq_file} 2> {log}"

rule fastqc_trimmed:
    input:
        fastq_file="output/trimmed/{sample}_{read}_trimmed.fastq.gz"
    output:
        zip_files="output/fastqc/{sample}_{read}_trimmed_fastqc.zip",
        html_files="output/fastqc/{sample}_{read}_trimmed_fastqc.html"
    threads: lambda wildcards, attempt: min(attempt * 2, 8)
    log:
        "output/fastqc/{sample}_{read}_trimmed_fastqc.log"
    container:
        "docker://staphb/fastqc"
    shell:
        "fastqc -t {threads} --outdir output/fastqc {input.fastq_file} 2> {log}"



    