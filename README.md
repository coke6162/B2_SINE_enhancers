# Co-option of B2 SINE elements as interferon-inducible enhancers in mouse

## Data availability:
All raw and processed sequencing data generated in this study have been submitted to the NCBI Gene Expression Omnibus (GEO) with accession number [GSE202574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202574).

## Publicly available data used:
List of publicly available data used in this study:
* [Piccolo, V., Curina, A., Genua, M. et al. Opposing macrophage polarization programs show extensive epigenomic and transcriptional cross-talk. Nat Immunol 18(5), 530–540 (2017).](https://doi.org/10.1038/ni.3710) RNA-seq & ChIP-seq (GSE84520)
* [Platanitis, E., Demiroz, D., Schneller, A. et al. A molecular switch from STAT2-IRF9 to ISGF3 underlies interferon-induced gene transcription. Nat Commun 10(1), 2921 (2019).](https://doi.org/10.1038/s41590-018-0184-1) RNA-seq & ChIP-seq (GSE115435)
* [Cuartero, S., Weiss, F.D., Dharmalingam, G. et al. Control of inducible gene expression links cohesin to hematopoietic progenitor self-renewal and differentiation. Nat Immunol 19(9), 932–941 (2018).](https://www.nature.com/articles/s41590-018-0184-1) RAD21 ChIP-seq (SRR6492207)
* [Gualdrini, F., Polletti, S., Simonatto, M., et al. H3K9 trimethylation in active chromatin restricts the usage of functional CTCF sites in SINE B2 repeats. Genes Dev 36(7-8), 414-432 (2022).](https://doi.org/10.1101%2Fgad.349282.121) CTCF ChIP-seq (SRR17090500, SRR17090494)
* [Platanitis, E., Gruener, S., Ravi Sundar Jose Geetha, A., et al. Interferons reshape the 3D conformation and accessibility of macrophage chromatin. iScience 25(3), (2022).](https://doi.org/10.1016/j.isci.2022.103840) ATAC-seq and Hi-C (PRJNA694816)
* [Qiao, Y., Kang, K., Giannopoulou, E., et al. IFN-γ Induces Histone 3 Lysine 27 Trimethylation at a Small Subset of Promoters to Stably Silence Gene Expression in Human Macrophages. Cell Rep. 16(12), 3121-3129 (2016).](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5079287/) Human RNA-seq (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84692)

## UCSC Genome Browser Session:
bigWig files for all samples analyzed in this study may be visualized on the UCSC Genome Browser [here](https://genome.ucsc.edu/s/coke6162/B2_SINE_enhancers_Horton_et_al).

## Programs used:
List of programs used for all analyses:
* BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
* FastQC v0.11.8 (https://github.com/s-andrews/FastQC)
* MultiQC v1.7 (https://github.com/ewels/MultiQC)
* HISAT2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
* Samtools v1.14 (http://www.htslib.org/)
* Subread v1.6.2 (http://subread.sourceforge.net/)
* DESeq2 v1.26.0 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* TEtranscripts v2.1.4 (https://github.com/mhammell-laboratory/TEtranscripts)
* MACS2 v2.1.1 (https://pypi.org/project/MACS2/)
* BWA v0.7.15 (https://github.com/lh3/bwa)
* Picard v2.6.0 (https://broadinstitute.github.io/picard/)
* MEME Suite v5.4.1 (https://meme-suite.org/meme/)
* Singularity v3.1.1 (https://github.com/hpcng/singularity)
* deepTools v3.5.1 (https://deeptools.readthedocs.io/en/develop/index.html)
* bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
* GIGGLE v0.6.3 (https://github.com/ryanlayer/giggle)
* FIMO v5.4.1 (https://meme-suite.org/meme/)
* bedGraphToBigWig v4 (http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)
* Bowtie2 v2.2.9 (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* Cufflinks v2.2.1 (http://cole-trapnell-lab.github.io/cufflinks/)
* Stringtie v1.3.3b (https://ccb.jhu.edu/software/stringtie/)
* Salmon v1.9.0 (https://combine-lab.github.io/salmon/)
* Pairix v0.3.7 (https://github.com/4dn-dcic/pairix)
* Juicer v1.6 (https://github.com/aidenlab/juicer)
* Activity-by-Contact Model v0.2.2 (https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction)

## Regulatory activity of B2_Mm2 in innate immunity
ChIP-seq and RNA-seq data in murine primary bone marrow derived macrophages (BMDMs) were downloaded from publicly available datasets and processed as described below. All data were aligned to mm10. 

#### 1. Identify interferon-inducible genes and transposon families

RNA-seq reads were assigned to gene annotation using Gencode vM19, and interferon stimulated genes (ISGs) were identified using DESeq2 comparing BMDMs stimulated with interferon gamma (IFNG) relative to untreated. Family-level transposable element (TE) expression was determined by realigning RNA-seq reads using hisat2, allowing multimappers (see [hisat2_k100.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2_k100.sbatch)). Reads were assigned to TE families using TEtranscripts with a custom GTF annotation file derived from Dfam annotation (see [generate_TEtranscripts_gtf.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/generate_TEtranscripts_gtf.sbatch)). IFNG-inducible TE families were identified using DESeq2.

**Full BMDM RNA-seq Workflow:**
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/multiqc.sbatch)
4. [hisat2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2.sbatch)
5. [featureCounts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/featureCounts.sbatch)
6. [DESeq2_genes.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/DESeq2_genes.R)
7. [extract_top_750_ISGs.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/extract_top_750_ISGs.sh), [extract_top_750_IRGs.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/extract_top_750_IRGs.sh), [extract_random_750_nonresponsive_genes.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/extract_random_750_nonresponsive_genes.sh)
8. [Dicer1_expression_bar.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/Dicer1_expression_bar.R)

**For TEtranscripts, realign bams to allow multiple alignments per read:**
1. [hisat2_k100.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/hisat2_k100.sbatch)
2. [generate_TEtranscripts_gtf.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/generate_TEtranscripts_gtf.sbatch)
3. [TEtranscripts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/TEtranscripts.sbatch)
4. [DESeq2_TEtranscripts.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/DESeq2_TEtranscripts.R)
5. [TEtranscripts_B2_bar_expression.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_BMDM/TEtranscripts_B2_bar_expression.R)

#### 2. Identify STAT1-bound regions and test for family-level TE enrichment

GIGGLE was used to create a database of all TE families in the mm10 mouse genome using Dfam v2.0 annotation. Results were filtered according to the reported odds ratio across H3K27ac and STAT1 peak regions (see [filter_giggle_results.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/filter_giggle_results.sh)). Predicted IFNG enhancer-TE associations were plotted as a bubble plot (see [giggle_bubbles.py](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_bubbles.py)). 

We then identified a subset of B2 SINE elements that are bound by STAT1 in IFNG-stimulated BMDMs. Proximity to the nearest interferon stimulated gene (ISG), interferon repressed gene (IRG), and nonresponsive gene was determined for each STAT1-bound element and plotted as a histogram (see [https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_piccolo_IFNG_4h_vs_UT.R](nearest_neighbor_histograms_piccolo_IFNG_4h_vs_UT.R)). We additionally plotted H3K27ac, STAT1, CTCF, and RAD21 ChIP-seq signal as well as predicted binding sites over B2_Mm2 elements as a heatmap (see [B2_Mm2_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/B2_Mm2_heatmap.sbatch)).

**A typical ChIP-seq workflow to call peaks looks like this:**
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/multiqc.sbatch)
4. [bwa.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bwa.sbatch)
5. [remove_duplicates.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/remove_duplicates.sbatch)
6. [macs2_piccolo.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_piccolo.sbatch), [macs2_platanitis.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_platanitis.sbatch), [macs2_cuartero.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_cuartero.sbatch), [macs2_gualdrini.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/macs2_gualdrini.sbatch)
7. [bdg_to_bigwig.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bdg_to_bigwig.sbatch)
8. [intersect_peak_replicates.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/intersect_peak_replicates.sbatch)
9. [xstreme.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/xstreme.sbatch)

**Repeat enrichment analysis:**
1. [giggle_index.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_index.sbatch)
2. [giggle_search](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_search.sbatch)
3. [filter_giggle_results.sh](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/filter_giggle_results.sh)
4. [giggle_bubbles.py](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/giggle_bubbles.py)
5. [overlap_TEs_STAT1_summits.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/overlap_TEs_STAT1_summits.sbatch)

**Nearest neighbor analysis:**
1. [get_overlapping_B2_nearest_neighbor.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/get_overlapping_B2_nearest_neighbor.sbatch)
2. [bedtools_closest_nearest_neighbor.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/bedtools_closest_nearest_neighbor.sbatch)
3. [nearest_neighbor_histograms_piccolo_IFNG_4h_vs_UT.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_piccolo_IFNG_4h_vs_UT.R), [nearest_neighbor_histograms_piccolo_IFNG_2h_vs_UT.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_piccolo_IFNG_2h_vs_UT.R), [nearest_neighbor_histograms_platanitis_IFNG_2h_vs_UT.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/nearest_neighbor_histograms_platanitis_IFNG_2h_vs_UT.R)

**Visualize ChIP-seq signal over B2_Mm2 as a heatmap:**
1. [get_overlapping_B2_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/get_overlapping_B2_heatmap.sbatch)
2. [B2_Mm2_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/B2_Mm2_heatmap.sbatch)

**Assess distribution of p-values for predicted GAS motifs that overlap B2 SINE elements as a box-and-whisker plot:**
1. [get_overlapping_B2_box_whisker.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/get_overlapping_B2_box_whisker.sbatch)
2. [fimo_pval_1.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/fimo_pval_1.sbatch)
3. [get_fimo_pval_box_whisker.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/get_fimo_pval_box_whisker.sbatch)
4. [B2_GAS_whisker.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ChIPseq_BMDM/B2_GAS_whisker.R)

#### 3. Assess STAT1 and CTCF binding over B2 SINE subfamilies

We identified putative STAT1 and CTCF binding sites for the mm10 mouse genome assembly genome-wide using FIMO. B2 elements were "expanded" such that the coordinates are based on "full-length" boundaries relative to the consensus. Predicted motifs over all annotated B2 SINE elements were plotted as a heatmap (see [B2_Mm2_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm2_motif_heatmap.sbatch)). 

**Workflow:**
1. [fimo.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/fimo.sbatch)
2. [convert_fimo_txt_to_bw.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/convert_fimo_txt_to_bw.sbatch)
3. [run_createExpandedRepeatFile.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/run_createExpandedRepeatFile.sbatch)
4. [B2_Mm2_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm2_motif_heatmap.sbatch), [B2_Mm1a_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm1a_motif_heatmap.sbatch), [B2_Mm1t_motif_heatmap.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/motif_analysis/B2_Mm1t_motif_heatmap.sbatch)

#### 4. Predict enhancer-gene contacts using the Activity-by-Contact Model

We reanalyzed publicly available ATAC-seq and Hi-C data from IFNG-stimulated murine BMDMs and subsequently ran the Activity-by-Contact Model with H3K27ac ChIP-seq and RNA-seq generated from previous analyses to predict enhancer-gene contacts.

**ATAC-seq Workflow:**
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/multiqc.sbatch)
4. [bowtie2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/bowtie2.sbatch)
5. [remove_duplicates.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/remove_duplicates.sbatch)
6. [shift_fragments.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/shift_fragments.sbatch)
7. [macs2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ATACseq_BMDM/macs2.sbatch)

**Hi-C Workflow:**
1. [bwa.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/bwa.sbatch)
2. [pairtools_parse.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_parse.sbatch)
3. [pairtools_sort.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_sort.sbatch)
4. [pairtools_merge.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_merge.sbatch)
5. [pairtools_dedup.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_dedup.sbatch)
6. [pairtools_filter.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/pairtools_filter.sbatch)
7. [run_generate_site_positions.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/run_generate_site_positions.sbatch)
8. [run_addfrag2pairs.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/run_addfrag2pairs.sbatch)
9. [juicer_pre.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/juicer_pre.sbatch)
10. [juicebox_dump.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/HiC_BMDM/juicebox_dump.sbatch)

**ABC Workflow:**
1. [call_candidate_regions.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/call_candidate_regions.sbatch)
2. [collapse_gene_boundaries.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/collapse_gene_boundaries.sbatch)
3. [find_neighborhoods.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/find_neighborhoods.sbatch)
4. [predict_enhancers.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/predict_enhancers.sbatch)
5. [subset_and_intersect_enhancers.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/ABC_BMDM/subset_and_intersect_enhancers.sbatch)

#### 5. Comparative analysis assessing regulatory contibutions of B2_Mm2 on innate immunity in mouse

We reanalyzed publicly available RNA-seq data from human CD14+ monocytes stimulated with IFNG for 24 hours to define a set of human ISGs, from which human-to-mouse one-to-one orthologs were identified. Mouse and human-to-mouse ISGs were binned according to species specificity, and the nearest STAT1-bound B2_Mm2 elements relative to each ISG or top putative B2_Mm2 enhancer predicted to interact with an ISG were identified.

**Human RNA-seq Analysis Workflow:**
1. [bbduk.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/bbduk.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/multiqc.sbatch)
4. [hisat2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/hisat2.sbatch)
5. [featureCounts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/featureCounts.sbatch)
6. [DESeq2.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/DESeq2.R)

**Orthology Analysis Workflow:**
1. [generate_binary_matrix.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/generate_binary_matrix.sbatch)
2. [identify_nearest_STAT1_B2_Mm2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/identify_nearest_STAT1_B2_Mm2.sbatch)
3. [identify_interacting_ABC_B2_Mm2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/orthology_analysis/identify_interacting_ABC_B2_Mm2.sbatch)

#### 6. CRISPR-mediated deletion of B2_Mm2.Dicer1

We generated J774A.1 clones harboring a deletion for a B2_Mm2 element intronic to the *Dicer1* gene. Changes in gene expression were quantified using qPCR (see [qPCR_bargraph.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/qPCR_bargraph.R)) and RNA-seq. 

**Mutant J774A.1 RNA-seq Workflow:**
1. [bbduk_PE.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/bbduk_PE.sbatch)
2. [fastqcreport.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/fastqcreport.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/multiqc.sbatch)
4. [hisat2_PE.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/hisat2_PE.sbatch)
5. [merge_bams.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/merge_bams.sbatch)
6. [bamTobw.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/bamTobw.sbatch)
7. [featureCounts.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/fastqcreport.sbatch)
8. [deseq2.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/deseq2.R)
9. [distance_plots.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/distance_plots.R)
10. [normalized_CPM_plots.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/RNAseq_J774/normalized_CPM_plots.R)

To assess whether deletion of the B2_Mm2 element changes relative isoform abundances, we ran Stringtie to assemble novel transcripts and subsequently performed differential expression analysis at the transcript level.

**Isoform Expression Workflow:**
1. [stringtie_assemble.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/stringtie_assemble.sbatch)
2. [stringtie_merge.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/stringtie_merge.sbatch), [gtf_list.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/gtf_list.txt)
3. [salmon_index_decoy.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/salmon_index_decoy.sbatch), [decoys.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/decoys.txt)
4. [salmon_quant.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/salmon_quant.sbatch)
5. [DESeq2.R](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/DESeq2.R), [sample_names.txt](https://github.com/coke6162/B2_SINE_enhancers/blob/main/isoform_expression_J774/sample_names.txt)

We additionally performed CUT&TAG (H3K27ac, STAT1, POLR2A) on wild-type J774A.1 cells and J774A.1 cells harboring a deletion for B2_Mm2.Dicer1. 

**CUT&TAG Workflow:**
1. [bbduk_PE.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/bbduk_PE.sbatch)
2. [fastqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/fastqc.sbatch)
3. [multiqc.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/multiqc.sbatch)
4. [bwa_batch.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/bwa_batch.sbatch)
5. [MACS2.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/MACS2.sbatch)
6. [bdg_bw.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/bdg_bw.sbatch)
7. [calculate_FRIP_score.sbatch](https://github.com/coke6162/B2_SINE_enhancers/blob/main/CUTnTAG_J774/calculate_FRIP_score.sbatch)