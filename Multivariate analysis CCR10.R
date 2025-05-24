#––– libraries
library(tidyverse)
library(DESeq2)
library(gprofiler2)
library(survival)
library(survminer)
library(forestmodel)

#––– 1) read count matrix
raw_counts <- read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
counts_mat <- as.matrix(raw_counts %>% select(-Gene))
rownames(counts_mat) <- raw_counts$Gene

#––– 2) read and build translocation + CNA metadata (include del1p22!)
trans_df <- read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_tumor_only_mm_igtx_pairoscope.tsv")
cna_df   <- read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_gatk_cna_seqFISH.tsv")

trans_meta <- trans_df %>%
  transmute(
    sample_id = SAMPLE,
    patient   = str_extract(SAMPLE, "[0-9]{4}") %>% paste0("MMRF_", .),
    visit     = str_extract(SAMPLE, "_[0-9]_") %>% str_extract("[0-9]") %>% as.integer(),
    type      = str_extract(SAMPLE, "_[A-Z]{2}_") %>% str_extract("[A-Z]+"),
    relapsed  = if_else(visit == 1, 0L, 1L),
    NDS2      = NSD2_CALL,
    CCND1     = CCND1_CALL
  )

cna_meta <- cna_df %>%
  transmute(
    sample_id = SAMPLE,
    del17p13  = SeqWGS_Cp_17p13_20percent,
    Gain1q21  = SeqWGS_Cp_1q21_20percent,
    del1p22   = SeqWGS_Cp_1p22_20percent          # ← NEW
  )

cyto_combined <- trans_meta %>%
  filter(type == "BM") %>%
  left_join(cna_meta, by = "sample_id") %>%
  mutate_at(vars(relapsed, NDS2, CCND1, del17p13, Gain1q21, del1p22),
            ~ factor(.x, levels = c(0, 1)))

#––– 3) subset your count matrix to those same BM samples
common   <- intersect(colnames(counts_mat), cyto_combined$sample_id)
counts_mat    <- counts_mat[, common]
cyto_combined <- cyto_combined %>%
  filter(sample_id %in% common) %>%
  column_to_rownames("sample_id")

#––– 4) DESeq2 + VST
dds    <- DESeqDataSetFromMatrix(countData = counts_mat,
                                 colData   = cyto_combined,
                                 design    = ~1)
dds    <- estimateSizeFactors(dds)
vsd    <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

#––– 5) map ENSG → HGNC
mapping <- gconvert(query      = rownames(vst_mat),
                    organism   = "hsapiens",
                    target     = "HGNC",
                    mthreshold = 1,
                    filter_na  = FALSE)

vsd_df <- as.data.frame(vst_mat) %>%
  rownames_to_column("ensembl_id") %>%
  mutate(gene_symbol = mapping$target[match(ensembl_id, mapping$input)])

#––– 6) pull out CCR10
if (!"CCR10" %in% vsd_df$gene_symbol) {
  stop("CCR10 not found in your VST matrix via the gconvert mapping.")
}
CCR10_df <- vsd_df %>%
  filter(gene_symbol == "CCR10") %>%
  select(-ensembl_id, -gene_symbol) %>%
  pivot_longer(cols      = everything(),
               names_to  = "sample_id",
               values_to = "CCR10_expr")

#––– 7) read & clean survival
surv <- read_tsv("MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL.tsv") %>%
  transmute(
    patient = PUBLIC_ID,
    status  = if_else(!is.na(deathdy) & deathdy > 1, 1, 0),
    time    = coalesce(deathdy, lstalive, lvisitdy)
  )

#––– 8) build the cox_df including del1p22
cox_df <- CCR10_df %>%
  left_join(cyto_combined %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  left_join(surv, by = "patient") %>%
  select(time, status, CCR10_expr, relapsed, NDS2, CCND1,
         Gain1q21, del17p13, del1p22) %>%         # ← now includes del1p22
  drop_na()

#––– 9) fit & summarize multivariate Cox (with del1p22)
cox_mod <- coxph(
  Surv(time, status) ~ CCR10_expr +
    relapsed + NDS2 + CCND1 + Gain1q21 + del17p13 + del1p22,
  data = cox_df
)
summary(cox_mod)

#––– 10) forest plot
forest_model(cox_mod) +
  ggtitle("Multivariate Cox (now includes del1p22)")

#––– 11) example: CCR10 × del1p22 interaction
cox_int1p <- coxph(
  Surv(time, status) ~ CCR10_expr * del1p22 +
    relapsed + NDS2 + CCND1 + Gain1q21 + del17p13,
  data = cox_df
)
anova(cox_mod, cox_int1p, test = "Chisq")
summary(cox_int1p)

#––– 12) visualize that interaction
med_CCR10 <- median(cox_df$CCR10_expr, na.rm=TRUE)
newdata_1p <- expand.grid(
  CCR10_expr = c(med_CCR10 - 1, med_CCR10 + 1),
  del1p22   = levels(cox_df$del1p22),
  relapsed  = levels(cox_df$relapsed)[1],
  NDS2      = levels(cox_df$NDS2)[1],
  CCND1     = levels(cox_df$CCND1)[1],
  Gain1q21  = levels(cox_df$Gain1q21)[1],
  del17p13  = levels(cox_df$del17p13)[1]
) %>%
  mutate(group = paste0(ifelse(CCR10_expr<med_CCR10,"CCR10 low,","CCR10 high,"),
                        ifelse(del1p22=="1"," del1p22+"," del1p22-")))

fit_1p <- survfit(cox_int1p, newdata=newdata_1p)
ggsurvplot(fit_1p,
           data        = newdata_1p,
           legend.labs = newdata_1p$group,
           palette     = RColorBrewer::brewer.pal(4,"Dark2"),
           pval        = TRUE,
           conf.int    = FALSE,
           xlab        = "Days",
           ylab        = "Survival Probability",
           title       = "Interaction: CCR10 × del1p22")

