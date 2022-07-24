All bash scripts were written to be run using SLURM on the HPC cluster at the University of Colorado Boulder. Some scripts are written to be run in parallel as job arrays.

Accession numbers and corresponding sample names are provided in [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/sample_names.txt). These scripts assume that all samples are named as described in [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/sample_names.txt). Samples may be download from SRA using the provided SRR IDs using sra-tools (script not provided). For simplicity all input and output files are written to the same directory.

A typical ChIP-seq workflow to call peaks looks like this:
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/multiqc.sbatch)
4. [bwa.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bwa.sbatch)
5. [remove_duplicates.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/remove_duplicates.sbatch)
6. [macs2_piccolo.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_piccolo.sbatch), [macs2_platanitis.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_platanitis.sbatch), [macs2_cuartero.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_cuartero.sbatch), [macs2_gualdrini.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_gualdrini.sbatch)
7. [intersect_peak_replicates.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/intersect_peak_replicates.sbatch)
8. [xstreme.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/xstreme.sbatch)

Repeat enrichment analysis:
1. [giggle_index.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_index.sbatch)
2. [giggle_search](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_search.sbatch)
3. [filter_giggle_results.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/filter_giggle_results.sh)
4. [giggle_bubbles.py](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_bubbles.py)

Nearest neighbor analysis:
1. [get_overlapping_B2_nearest_neighbor.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/get_overlapping_B2_nearest_neighbor.sbatch)
2. [bedtools_closest.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bedtools_closest.sbatch)
3. [nearest_neighbor_histograms_piccolo_IFNG_4h_vs_UT.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_piccolo_IFNG_4h_vs_UT.R), [nearest_neighbor_histograms_piccolo_IFNG_2h_vs_UT.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_piccolo_IFNG_2h_vs_UT.R), [nearest_neighbor_histograms_platanitis_IFNG_2h_vs_UT.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_platanitis_IFNG_2h_vs_UT.R)

deepTools heatmap:
1. get overlapping B2_Mm2 by STAT1/CTCF occupancy
2. deepTools heatmap

Need following B2 overlap files:
* deepTools heatmap - B2_Mm2 bound only by STAT1, B2_Mm2 bound only by CTCF, B2_Mm2 bound by CTCF and STAT1, randomly shuffled unbound B2_Mm2

Bonus files:
* Overlap H3K27ac & STAT1 (referenced in text)
* Overlap all TEs with STAT1 summits & motifs (supp table S3)
* Nearest neighbor repeat, using log2FC cutoff of > 0.5 (supp table S5)
* Box-and-whisker (supp figure S5)
* Differential AME analysis (supp table S6) - to be removed?

Required packages:
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* BWA v0.7.15 (https://github.com/lh3/bwa)
* Samtools v1.10 (http://www.htslib.org/)
* Picard v2.6.0 (https://broadinstitute.github.io/picard/)
* MACS v2.1.1 (https://pypi.org/project/MACS2/)
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
* XSTREME v5.4.1 (https://meme-suite.org/meme/)
* GIGGLE v0.6.3 (https://github.com/ryanlayer/giggle)
* deepTools v3.5.1 (https://deeptools.readthedocs.io/en/develop/index.html)

Note that the soft masked mm10 assembly file may be downloaded via:
* mm10.fa - http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz