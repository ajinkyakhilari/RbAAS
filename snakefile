import os

SAMPLES_DIR = config["sample_dir"]
OUTPUT_DIR = config["output_dir"]
BARCODES = [d for d in os.listdir(SAMPLES_DIR) if os.path.isdir(os.path.join(SAMPLES_DIR, d))]

rule all:
  input:
    expand("{output_dir}/{barcode}/{barcode}.final.consensus.fasta", barcode=BARCODES, output_dir=OUTPUT_DIR)

rule merge_fastq:
  input:
    dummy=lambda wildcards: os.path.join(config["sample_dir"], wildcards.barcode)
  output:
    merged="{output_dir}/{barcode}/{barcode}.fastq"
  conda:
    "envs/fastp.yaml"
  shell:
    """
    mkdir -p {wildcards.output_dir}/{wildcards.barcode}
    cat {input}/*.fastq > {output.merged}
    """

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
    minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} |
      samtools view -bS -F 4 - | 
      samtools sort -@ {threads} -o {output.bam}
    samtools index {output.bam}
    """

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
    ivar trim -i {input.bam} -b {input.primers} -e -m 1 -s 4 -q 0 -p {wildcards.output_dir}/{wildcards.barcode}/{wildcards.barcode}.primerclipped
    samtools sort -@ {threads} {wildcards.output_dir}/{wildcards.barcode}/{wildcards.barcode}.primerclipped.bam \
      -o {output.trimmed}
    samtools index {output.trimmed}
    """

rule clair3_call:
  input:
    bam="{output_dir}/{barcode}/{barcode}.primerclipped.sorted.bam",
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
                  --output {wildcards.output_dir}/{wildcards.barcode}/clair3_output \
                  --include_all_ctgs --haploid_sensitive
    tabix -p vcf {output.vcf} -f
    """

rule mask_variants:
  input:
    vcf="{output_dir}/{barcode}/clair3_output/merge_output.vcf.gz",
    bam="{output_dir}/{barcode}/{barcode}.primerclipped.sorted.bam"
  output:
    mask="{output_dir}/{barcode}/{barcode}.mask.bed"
  conda:
    "envs/bcftools.yaml"
  shell:
    """
    bcftools query -f '%CHROM\t%POS0\t%POS\n' \
      -i 'QUAL<10 || FORMAT/GQ<10 || FORMAT/AF<0.5' {input.vcf} > \
      {wildcards.output_dir}/{wildcards.barcode}/{wildcards.barcode}.variant_mask.bed

    samtools depth -a {input.bam} | awk '$3 < 20 {{print $1, $2-1, $2}}' OFS="\t" > \
      {wildcards.output_dir}/{wildcards.barcode}/{wildcards.barcode}.coverage_mask.bed

    cat {wildcards.output_dir}/{wildcards.barcode}/{wildcards.barcode}.variant_mask.bed \
        {wildcards.output_dir}/{wildcards.barcode}/{wildcards.barcode}.coverage_mask.bed | \
        sort -k1,1 -k2,2n | uniq > {output.mask}
    """

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
    bcftools consensus -f {input.ref} -m {input.mask} -o {output.consensus} {input.vcf}
    """

