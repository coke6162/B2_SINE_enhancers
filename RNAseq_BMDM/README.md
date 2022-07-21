All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

Accession numbers and corresponding sample names are provided in sample_names.txt. These scripts assume that all samples are named as described in sample_names.txt. For simplicity all input and output files are written to the same directory.

A typical RNA-seq workflow to call differentially expressed genes looks like this:
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/multiqc.sbatch)
4. [hisat2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2.sbatch)
5. [featureCounts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/featureCounts.sbatch) (references bam_order.txt)
6. [DESeq2_genes.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/DESeq2_genes.R) (references gencode_vM18_crossref.txt and gencode_vM18_tss.bed)
7. Get top 750 ISGs, IRGs, nonresponsive genes
8. Dicer1 gene expression count plots

A typical RNA-seq workflow to call differentially expressed TE families looks like this:
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/multiqc.sbatch)
4. [hisat2_k100.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2_k100.sbatch)
5. [generate_TEtranscripts_gtf.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/generate_TEtranscripts_gtf.sbatch)
6. [TEtranscripts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/TEtranscripts.sbatch)
7. [DESeq2_TEtranscripts.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/DESeq2_TEtranscripts.R)
8. B2_Mm2 expression count plots]()

Required packages:
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* HISAT2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
* Samtools v1.10 (http://www.htslib.org/)
* Subread v1.6.2 (http://subread.sourceforge.net/)
* DESeq2 v1.26.0 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* TEtranscripts v2.1.4 (https://github.com/mhammell-laboratory/TEtranscripts)

For a list of all R packages used in this analysis, see session_info.txt.

Note that the soft masked mm10 assembly and Gencode vM18 gene annotation file may be downloaded via:
* mm10.fa - http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
* gencode.vM19.annotation.gtf - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.annotation.gtf.gz