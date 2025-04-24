configfile: "config.yaml"

import os

# Directories and barcodes
SAMPLES_DIR = config["sample_dir"]
OUTPUT_DIR = config["output_dir"]
BARCODES = [d for d in os.listdir(SAMPLES_DIR) if os.path.isdir(os.path.join(SAMPLES_DIR, d))]

# Flag to enable or disable primer trimming
TRIM_PRIMERS = config.get("trim_primers", False)

# Default target: all final consensus FASTA files
rule all:
    input:
        expand(
            "{output_dir}/{barcode}/{barcode}.final.consensus.fasta",
            barcode=BARCODES,
            output_dir=OUTPUT_DIR
        )

# 1. Merge all FASTQ files per barcode
rule merge_fastq:
    input:
        dummy=lambda wc: os.path.join(SAMPLES_DIR, wc.barcode)
    output:
        merged="{output_dir}/{barcode}/{barcode}.fastq"
    conda:
        "envs/fastp.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/{wildcards.barcode}
        cat {input}/*.fastq > {output.merged}
        """

# 2. Quality filtering with fastp
rule quality_filter:
    input:
        "{output_dir}/{barcode}/{barcode}.fastq"
    output:
        filt="{output_dir}/{barcode}/{barcode}_filt.fastq"
    params:
        quality=config["quality"],
        min_length=config["min_length"],
        max_length=config["max_length"]
    conda:
        "envs/fastp.yaml"
    threads: config["threads"]
    shell:
        """
        fastp -i {input} -o {output.filt} \
              -q {params.quality} -l {params.min_length} --max_len1 {params.max_length} \
              --thread {threads}
        """

# 3. Alignment to reference with minimap2 & indexing
rule align_minimap2:
    input:
        fastq="{output_dir}/{barcode}/{barcode}_filt.fastq",
        ref=config["ref"]
    output:
        bam="{output_dir}/{barcode}/{barcode}.sorted.bam"
    conda:
        "envs/minimap2.yaml"
    threads: config["threads"]
    shell:
        """
        samtools faidx {input.ref}
        minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} | \
          samtools view -bS -F 4 - | \
          samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# 4. Optional primer trimming with ivar
if TRIM_PRIMERS:
    rule trim_primers:
        input:
            bam="{output_dir}/{barcode}/{barcode}.sorted.bam",
            primers=config["primer_bed"]
        output:
            trimmed="{output_dir}/{barcode}/{barcode}.primerclipped.sorted.bam"
        conda:
            "envs/ivar.yaml"
        threads: config["threads"]
        shell:
            """
            ivar trim -i {input.bam} -b {input.primers} -e -m 1 -s 4 -q 0 \
              -p {OUTPUT_DIR}/{wildcards.barcode}/{wildcards.barcode}.primerclipped
            samtools sort -@ {threads} \
              {OUTPUT_DIR}/{wildcards.barcode}/{wildcards.barcode}.primerclipped.bam \
              -o {output.trimmed}
            samtools index {output.trimmed}
            """

# Helper function to select BAM suffix based on trimming

def bam_input(wc):
    suffix = ".primerclipped.sorted.bam" if TRIM_PRIMERS else ".sorted.bam"
    return os.path.join(OUTPUT_DIR, wc.barcode, wc.barcode + suffix)

# 5. Variant calling with Clair3
rule clair3_call:
    input:
        bam=bam_input,
        ref=config["ref"]
    output:
        vcf="{output_dir}/{barcode}/clair3_output/merge_output.vcf.gz"
    conda:
        "envs/clair3.yaml"
    threads: config["threads"]
    shell:
        """
        MODEL_PATH=models/r941_prom_sup_g5014
        run_clair3.sh --bam_fn {input.bam} \
                      --ref_fn {input.ref} \
                      --model_path $MODEL_PATH \
                      --threads {threads} \
                      --platform ont \
                      --output {OUTPUT_DIR}/{wildcards.barcode}/clair3_output \
                      --include_all_ctgs --haploid_sensitive
        tabix -p vcf {output.vcf} -f
        """

# 6. Generate mask for low-quality and low-coverage regions
rule mask_variants:
    input:
        vcf="{output_dir}/{barcode}/clair3_output/merge_output.vcf.gz",
        bam=bam_input
    output:
        mask="{output_dir}/{barcode}/{barcode}.mask.bed"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools query -f '%CHROM\t%POS0\t%POS\n' \
          -i 'QUAL<10 || FORMAT/GQ<10 || FORMAT/AF<0.5' {input.vcf} > \
          {OUTPUT_DIR}/{wildcards.barcode}/{wildcards.barcode}.variant_mask.bed

        samtools depth -a {input.bam} | awk '$3 < 20 {{print $1, $2-1, $2}}' OFS="\t" > \
          {OUTPUT_DIR}/{wildcards.barcode}/{wildcards.barcode}.coverage_mask.bed

        cat {OUTPUT_DIR}/{wildcards.barcode}/{wildcards.barcode}.variant_mask.bed \
            {OUTPUT_DIR}/{wildcards.barcode}/{wildcards.barcode}.coverage_mask.bed | \
            sort -k1,1 -k2,2n | uniq > {output.mask}
        """

# 7. Build consensus sequence
rule generate_consensus:
    input:
        vcf="{output_dir}/{barcode}/clair3_output/merge_output.vcf.gz",
        mask="{output_dir}/{barcode}/{barcode}.mask.bed",
        ref=config["ref"]
    output:
        consensus="{output_dir}/{barcode}/{barcode}.final.consensus.fasta"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools consensus -f {input.ref} -m {input.mask} \
          -o {output.consensus} {input.vcf}
        """

# Note: Add `trim_primers: true` or `false` in config.yaml to toggle primer removal.

