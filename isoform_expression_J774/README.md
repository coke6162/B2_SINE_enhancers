All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

These scripts require the following files as input:
* Trimmed reads from the [J774 WT & B2_Mm2.Dicer1 KO RNA-seq analysis](https://github.com/coke6162/B2_SINE_enhancers/tree/main/RNAseq_J774)
* Aligned fragments from the [J774 WT & B2_Mm2.Dicer1 KO RNA-seq analysis](https://github.com/coke6162/B2_SINE_enhancers/tree/main/RNAseq_J774)

Workflow:
1. [stringtie_assemble.sbatch]()
2. [stringtie_merge.sbatch](), [gtf_list.txt]()
3. [salmon_index_decoy.sbatch](), [decoys.txt]
4. [salmon_quant.sbatch]()
5. [DESeq2.R](), [sample_names.txt]()

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
* gencode.vM18.annotation.gtf - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.annotation.gtf.gz