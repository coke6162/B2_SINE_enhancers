# B2_SINE_enhancers
WORK IN PROGRESS

## Scripts and files used in study:


## Data availability:
All raw and processed sequencing data generated in this study have been submitted to the NCBI Gene Expression Omnibus (GEO) with accession number [GSE202574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202574).

## Publicly available data used:
List of publicly available data used in this study.
* [Piccolo, V., Curina, A., Genua, M. et al. Opposing macrophage polarization programs show extensive epigenomic and transcriptional cross-talk. Nat Immunol 18, 530–540 (2017).](https://www.nature.com/articles/ni.3710) RNA-seq & ChIP-seq (GSE84520)
* [Platanitis, E., Demiroz, D., Schneller, A. et al. A molecular switch from STAT2-IRF9 to ISGF3 underlies interferon-induced gene transcription. Nat Commun 10, 2921 (2019).](https://www.nature.com/articles/s41467-019-10970-y) RNA-seq & ChIP-seq (GSE115435)
* [Cuartero, S., Weiss, F.D., Dharmalingam, G. et al. Control of inducible gene expression links cohesin to hematopoietic progenitor self-renewal and differentiation. Nat Immunol 19, 932–941 (2018).](https://www.nature.com/articles/s41590-018-0184-1) RAD21 ChIP-seq (SRR6492207)
* [Gualdrini, F., Polletti, S., Simonatto, M., et al. H3K9 trimethylation in active chromatin restricts the usage of functional CTCF sites in SINE B2 repeats. Genes Dev 36, 414-432 (2022)](http://genesdev.cshlp.org/content/early/2022/03/30/gad.349282.121) CTCF ChIP-seq (SRR17090500, SRR17090494)

## Programs used:
List of programs used for all analyses.
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* HISAT2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
* Samtools v1.10 (http://www.htslib.org/)
* Subread v1.6.2 (http://subread.sourceforge.net/)
* DESeq2 v1.26.0 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* TEtranscripts v2.1.4 (https://github.com/mhammell-laboratory/TEtranscripts)
* MACS2 v2.1.1 (https://pypi.org/project/MACS2/)
* BWA v0.7.15 (https://github.com/lh3/bwa)
* MEME Suite v5.4.1 (https://meme-suite.org/meme/)
* Singularity v3.1.1 (https://github.com/hpcng/singularity)
* deepTools v3.0.1 (https://deeptools.readthedocs.io/en/develop/index.html)
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)

## Regulatory activity of B2_Mm2 in innate immunity
ChIP-seq and RNA-seq data in murine primary bone marrow derived macrophages (BMDMs) were downloaded from publicly available datasets and processed as described below. All data were aligned to mm10. 

1. Identify interferon-inducible genes and transposon families

RNA-seq reads were assigned to gene annotation using Gencode vM19, and interferon stimulated genes (ISGs) were identified using DESeq2 comparing BMDMs stimulated with interferon gamma (IFNG) relative to untreated. Family-level transposable element (TE) expression was determined by realigning RNA-seq reads using hisat2, allowing multimappers (see [hisat2_k100.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2_k100.sbatch)). Reads were assigned to TE families using TEtranscripts with a custom GTF annotation file derived from Dfam annotation (see [generate_TEtranscripts_gtf.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/generate_TEtranscripts_gtf.sbatch)). IFNG-inducible TE families were identified using DESeq2.

2. Identify STAT1-bound regions and test for family-level TE enrichment


3. Assess STAT1 and CTCF binding over B2 SINE subfamilies

RNA-seq Workflow:
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/multiqc.sbatch)
4. [hisat2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2.sbatch)
5. [featureCounts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/featureCounts.sbatch) (references bam_order.txt)
6. [DESeq2_genes.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/DESeq2_genes.R) (references gencode_vM18_crossref.txt and gencode_vM18_tss.bed)
7. [extract_top_750_ISGs.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/extract_top_750_ISGs.sh), [extract_top_750_IRGs.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/extract_top_750_IRGs.sh), [extract_random_750_nonresponsive_genes.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/extract_random_750_nonresponsive_genes.sh)
8. [Dicer1_expression_bar.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/Dicer1_expression_bar.R)

For TEtranscripts, realign bams to allow multiple alignments per read:
1. [hisat2_k100.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2_k100.sbatch)
2. [generate_TEtranscripts_gtf.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/generate_TEtranscripts_gtf.sbatch)
3. [TEtranscripts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/TEtranscripts.sbatch)
4. [DESeq2_TEtranscripts.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/DESeq2_TEtranscripts.R)
5. [TEtranscripts_B2_bar_expression.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/TEtranscripts_B2_bar_expression.R)

ChIP-seq Workflow:

## B2_Mm2 transcription factor profiling
Details details.

Workflow:

## CRISPR-mediated deletion of B2_Mm2.Dicer1
Details details.

RNA-seq Workflow:

CUT&TAG Workflow: