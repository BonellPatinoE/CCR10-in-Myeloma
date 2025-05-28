#CCR10 in multiple myeloma

This repository contains two different single-cell RNA-sequencing analyses about CCR10 expression in primary samples from patients with multiple myeloma. CCR10 is a a chemokine receptor discovered at Wiita lab using cell surface capture proteomic approach as a novel therapeutic target in multiple myeloma (Nat Commun 13, 4121 2022). 

The first analysis was made from a public data previously published (GSE223060, Cancer Res. 2023 Apr 15;83(8):1214–33) and was conducted in R (v4.1.2) with Seurat (v4.3.0)2. Samples with more than 25% of plasma cells were selected for the analysis. (MMRF_1325,MMRF1537,   MMRF1640,   MMRF_2038,   MMRF_1267,   MMRF_1720,   MMRF_1505,   MMRF_2251,MMRF_2259).

Second analysis was conducted from a public available scRNAseq collection of scRNA-seq and scTCR-seq data from more than 200 samples from 182 patients with multiple myeloma, MGUS and SMM, and non-cancer controls (https://zenodo.org/records/13646014, Foster KA, et al http://medrxiv.org/lookup/doi/10.1101/2024.06.22.24309250). This data was primarily collected from  omnibus (GEO) Maura et al. under accession GSE161195 (Nat Cancer. 2023 Dec 1;4(12):1660–74, Nat Med. 2021 Mar 1;27(3):491–503), Bailur et al. (GSE163278, JCI Insight. 2019 Jun 6;4(11)), Oetjen et al. (GSE12022110, JCI Insight. 2018 Dec 6;3(23), Granja et al. (GSE139369, Nat Biotechnol. 2019 Dec 1;37(12):1458–65), Zavidij et al. (GSE124310, Nat Cancer. 2020 May 1;1(5):493–506), Kfoury   et   al.   (GSE143791, Cancer Cell. 2021 Nov 8;39(11):1464-1478.e8),   and   Zheng   et   al.   (GSE156728, Science (1979). 2021 Dec 17;374(6574)).   Using   a   paninmune   collection (panImmune.h5ad) it was analyzed it in Python using Scanpy (v1.9) through google colabs, using A100 GPU run time type. Each sample was categorized into “HD,” “MGUS,” “SMM,” or “MM” by matching sample
ID strings against pre‐defined regular expressions.

RNAseq data analysis about CCR10 expression in primary samples from coMMpass (MMRF) IA19A. This analysis was performed using DESeq2 (v1.34.0) package in R (v4.1.2) as well, using 'apeglm' for LFC shrinkage and “gprofiler” (v0.2.2) for for gene list functional enrichment analysis and namespace conversion. 

Bioinformatic analysis performed at Wiita Lab at the University of California, San Francisco (UCSF). 
