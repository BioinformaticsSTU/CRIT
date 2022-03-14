# CRIT

CRIT (**C**ircRNA **R**egulator **I**dentification **T**ool) is a pipeline based on a non-negative matrix factorization method to integrate various omics information and to identify regulating RBPs. 

# Table of contents

- [CRIT](#CRIT)
- [Table of contents](#Table-of-contents )
- [How CRIT works (in brief)](#How-CRIT-works-(in-brief) )
- [Before you begin](#Before-you-begin )
- [CRIT input files](#CRIT-input-files )
- [Running CRIT](#Running-CRIT )
- [CRIT output files](#CRIT-output-files )

# How CRIT works (in brief)

[(Back to top)](#table-of-contents)

CRIT adopts an NMF-based strategy consisting of a *"features collecting"* , a *"cluster identify"* and a *"Validate candidates"* component, based on the underlying assumptions that for a given gene to mediate a circRNA that gene must be: (1) Since the expression of regulators could affect the circRNAs expression, therefore regulators and some circRNA should have a strong expression correlation; (2) Regulators should have at least one of the following alterations, such as point mutation, copy number variation and abnormal mRNA expression in the cancer sample;(3) Above alterations on regulators have an impact on the survival time of patients; (4) Regulators' biological role should be closely related to the GO terms 'mRNA splicing, via spliceosome', which indicate the participation in the formation of circRNA; (5) Regulators can physically interact with some circRNAs.

**1.1. features collecting**. CRIT will collect features in five ways:

**ρm**: Based on expression data from MiOncoCirc compendium, for each gene, the mean value of absolute spearman's rank correlation coefficients of the top 50 gene-circRNAs pairs will be calculated, which captures the main influence of this gene on circRNAs expression.

**Zp**: From the TCGA project which consists of 19 batches of 12 cancer types, a Z score built from the dichotomized original data for each gene, which reflects the proportion of gene alterations in the samples of one cancer type or batch, which represents Z-score of gene alteration proportion.

**Zs**: Log-rank test between samples with or without gene alteration, which is implemented by the R package ‘survival’. The value for each gene was calculated from survival analysis, which reflects how significantly the alterations of the gene impact the patient survival data. Then each gene is assigned a Z score of survival analysis p-value.

**Zg**: Between GO term the gene annotated and with one functional term (GO:0000398, mRNA splicing, via spliceosome), is computed by R package 'GoSemSim', which is the Z-score of GO similarity scores between genes and selected GO term. 

**C value**: From the RNAct database, excluding the weak relationship with interaction score<30 of the protein-RNA interaction data, then align circRNA to RNA according to their chromosomal position annotations, which were extracted respectively from Gencode or CircBase. It is considered that a protein interacts with a circRNA if the position of the circRNA is completely aligned to the corresponding RNA.*C* value was defined as the number of circRNA with which each protein binds, and is used to indicate the strength of a protein interacting with circRNA. 

**1.2. cluster identify**. CRIT has been designed to identify circRNA regulators by using an NMF-based pipeline and pulling out genes that were collected in the previous research. We provide an .RData file containing a list of genes involved in 1344 RBPs. The user can of course prepare a similar .RData file containing genes relating to circRNA regulators or RBPs for input to CRIT.  OUr .RData file containing genes was prepared as follows:

- **(i)** Ensembl gene id: Ensembl gene IDs begin with ENS for Ensembl, and then a G for gene, such as ENSG00000148584.
- **(ii)** Hgnc symbol: The HGNC approves a gene name and symbol (short-form abbreviation) for each known human gene, such as A1CF.
- **(iii)** Entrez gene id: Identifier for a gene per the NCBI Entrez database, such as 29974.
- **(iv)** UniProt SwissProt: all genes annotated to the protein level, such as Q9NQ94.

Then the NMF-based method was applied for the features matrix.

**1.3. Validate candidates**. CRIT will also validate genes that appear in the cluster which was enriched in gold regulators. 
A list of gold regulators was collected from the previous research and assays based literature. Then multiple evidence from the Harmonizome database were collected to support candidate regulators‘ influence on circRNAs.

# Before you begin

[(Back to top)](#table-of-contents)

**2.1.** Download the necessary files to run CRIT and place them ALL locally into a directory of your choosing.   Your directory should contain the following files:

- **"GET_Prop_zscore.R"** (Z-Score of alteration in all samples)
- **"GET_Sur_zscore.R"** (Z-Score of each gene and the effect of the gene on the survival time of patients)
- **"GET_Gosim_zscore.R"** (Z-Score of similarity between GO term)
- **"GET_mRNAcirc_top50cor.R"** (mean value of top 50 pairs of circRNAs and genes)
- **"GET_C_value.R"** (numbers of circRNAs aligin to gene)
- **"GET_feauters_MATRIX.R"** (merge feauters to a matrix)
- **"NMF.R"** (NMF method to cluster candidate RBPs)
- **"RBP_analysis.R"** (candidate RBPs validation)
- **"RBP_list.RData"** (OPTIONAL: a list of 1344 candidate RBPs)

**2.2.** Finally, ensure you have downloaded the following R packages and dependencies:

- **clusterProfiler** (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
- **survival** (https://cran.r-project.org/web/packages/survival/survival.pdf)
- **org.Hs.eg.db** (https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
- **biomaRt** (https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- **data.table** (https://cran.r-project.org/web/packages/data.table/data.table.pdf)
- **GOSemSim** (http://bioconductor.org/packages/release/bioc/html/GOSemSim.html)
- **GOSim** (http://bioconductor.org/packages/release/bioc/html/GOSim.html)
- **readr** (https://cran.r-project.org/web/packages/readr/readr.pdf)
- **stringr** (https://cran.r-project.org/web/packages/stringr/stringr.pdf)
- **NMF** (https://cran.r-project.org/web/packages/NMF/NMF.pdf)
- **ggplot2** (https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf)
- **magrittr** (https://cran.r-project.org/web/packages/magrittr/magrittr.pdf)
- **dplyr** (https://cran.r-project.org/web/packages/dplyr/dplyr.pdf)
- **ComplexHeatmap** (https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt","survival","org.Hs.eg.db","clusterProfiler","GOSemSim","GOSim","readr","stringr","NMF","ggplot2","ComplexHeatmap")
install.packages("xlsx","data.table","dplyr","magrittr")
```

# CRIT input files

[(Back to top)](#table-of-contents)

The user must prepare following files to the upper directory named ./DATA/input folder which ./CODE folder already containing the GET_*.R and .RData files.  These input files must contain and be formatted as follows:

**3.1.** TCGA data(./DATA/input/TCGA/cancer/)

-  A tab-separated .txt file containing mRNA expression data of TCGA cancer types named **data_mRNA_median_Zscores.txt** or **data_RNA_Seq_v2_mRNA_median_Zscores.txt**.
-  A tab-separated .txt file containing CNV data of TCGA cancer types named **data_CNA.txt**.
-  A tab-separated .txt file containing mutation data of TCGA cancer types named **data_mutations_extended.txt.**

**3.2.** MiOncoCirc data(./DATA/input/MiOncoCirc)

The files below can be downloaded from https://mioncocirc.github.io/download/.

- **fpkm_matrix.csv** : gene expression file in FPKM 
- **v0.1.release.txt** : All circRNAs (filtered at 2 reads) in 2,000+ cancer samples in hg38
- **meta_update.xlsx** : All 2,000+ samples annotation

**3.3.** RNAct data(./DATA/input/RNAct)

- **catrapid_human_basic.zip** : interactome files contain human genome-wide catRAPID interaction prediction scores, which contains 20,778 canonical Human Reference Proteome proteins and 98,608 human GENCODE "basic" RNAs. This file can be download from https://rnact.crg.eu/download.
- **gencode.v37.annotation.gtf** : human RNA annotation file which can download from https://www.gencodegenes.org/human/release_37.html.
- **hsa_hg19_circRNA.txt** : human circRNA annotation file which can download from http://circbase.org/download/hsa_hg19_circRNA.txt.

**3.4.** Harmonizome data(./DATA/input/Harmonizome)

- **Physical interactions** (https://maayanlab.cloud/Harmonizome/dataset/Pathway+Commons+Protein-Protein+Interactions)
- **Cancer Gene Co-expression Modules** (https://maayanlab.cloud/Harmonizome/dataset/MSigDB+Cancer+Gene+Co-expression+Modules)
- **Curated Transcription Factor Targets** (https://maayanlab.cloud/Harmonizome/dataset/TRANSFAC+Curated+Transcription+Factor+Targets)

**3.5.** gold regulators list data(./DATA/input/)

- **golden_regulators.xlsx** : A .xlsx file contain gold regulators collected from previous researches.

# Running CRIT

[(Back to top)](#table-of-contents)

After downloading, preparing, and making local copies of the necessary files (see above), the user should specify their path parameters in the "**CRIT_settings.R**" file where indicated.  This file can then be saved and sourced, and CRIT will then run in its entirety.

# CRIT output files

[(Back to top)](#table-of-contents)

After successful completion, CRIT will write multiple data tables to file to output directory.  The output files are:

