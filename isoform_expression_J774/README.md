All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

These scripts require trimmed reads & aligned fragments from the [J774 WT & B2_Mm2.Dicer1 KO RNA-seq analysis](https://github.com/coke6162/B2_SINE_enhancers/tree/main/RNAseq_J774) as input.

Workflow:
1. [stringtie_assemble.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/stringtie_assemble.sbatch)
2. [stringtie_merge.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/stringtie_merge.sbatch), [gtf_list.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/gtf_list.txt)
3. [salmon_index_decoy.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/salmon_index_decoy.sbatch), [decoys.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/decoys.txt)
4. [salmon_quant.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/salmon_quant.sbatch)
5. [DESeq2.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/DESeq2.R), [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/sample_names.txt)

Required packages:
* Samtools v1.10 (http://www.htslib.org/)
* DESeq2 v1.26.0 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* Cufflinks v2.2.1 (http://cole-trapnell-lab.github.io/cufflinks/)
* Stringtie v1.3.3b (https://ccb.jhu.edu/software/stringtie/)
* Salmon v1.9.0 (https://combine-lab.github.io/salmon/)

Note that the soft masked mm10 assembly and Gencode vM18 gene annotation file may be downloaded via:
* [mm10.fa](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)
* [gencode.vM18.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.annotation.gtf.gz)