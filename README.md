# 🧬 nanoAmpRef: Reference-Based Amplicon Assembly Pipeline (Nanopore)

`nanoAmpRef` is a Snakemake-powered pipeline for **reference-based assembly of amplicons** sequenced using **Oxford Nanopore Technology (ONT)**. It performs quality filtering, alignment, primer trimming, variant calling, and masked consensus sequence generation — perfect for amplicon panels like *Plasmodium falciparum* drug resistance.

---

## 🔧 Key Features

- 🔬 Works with barcoded Nanopore amplicon reads
- 🧼 Quality filtering with `fastp`
- 🧭 Alignment to reference with `minimap2`
- ✂️ Primer trimming using `iVar`
- 🧬 Variant calling using `Clair3`
- 🩻 Low quality and low coverage region masking
- ✅ Final consensus FASTA file generation
- 📦 Fully reproducible with `--use-conda`

---

## 📁 Directory Structure

project/ ├── snakefile ├── envs/ │ ├── fastp.yaml │ ├── minimap2.yaml │ ├── ivar.yaml │ ├── clair3.yaml │ └── bcftools.yaml ├── run_pipeline.sh # optional wrapper └── fastq_pass/ ├── barcode01/ │ └── *.fastq └── barcode02/ └── *.fastq



---

## ⚙️ Installation

```bash
conda install -c bioconda -c conda-forge snakemake
```


## 🚀 Usage

## ✅ Basic Run

```bash
snakemake --cores 8 --use-conda \
  --config \
    sample_dir=/path/to/directory \
    output_dir=results \
    ref=/path/to/reference.fasta \
    primer_bed=path/to/primer.bed \
    quality=10 \
    min_length=150 \
    max_length=2000 \
    threads=8
```

## ✅ Run from another directory

```bash
snakemake --cores 8 --use-conda \
  --snakefile /full/path/to/snakefile \
  --directory /path/to/project \
  --config \
    sample_dir=/path/to/directory \
    output_dir=results \
    ref=/path/to/reference.fasta \
    primer_bed=path/to/primer.bed \
    quality=10 \
    min_length=150 \
    max_length=2000 \
    threads=8
```

## 🔬 Outputs per Barcode

    <output_dir>/<barcode>/<barcode>.fastq – merged input

    <barcode>.filt.fastq – quality-filtered

    <barcode>.primerclipped.sorted.bam – trimmed & aligned

    clair3_output/merge_output.vcf.gz – variant calls

    <barcode>.mask.bed – variant + coverage masking

    <barcode>.final.consensus.fasta – final consensus

## 📜 License

MIT License







