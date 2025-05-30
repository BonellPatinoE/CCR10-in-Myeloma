library(tidyverse)
library(tibble)
library(survminer)
library(survival)

################# Survival curve based on CCR10 expression (ENSG00000184451) ###################



## Read HtSeq data Survival based on 20% top and 20% bottom
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
#raw <- data.frame(Expression.Estimates...Gene.Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts)
ncol(raw)
raw<-raw%>%filter(Gene == "ENSG00000184451")
raw
counts <- as.matrix(raw[-1])
counts1<-data.frame(t(counts))
counts1
ncol(counts)



counts2<- tibble::rownames_to_column(counts1, "Patient")
colnames(counts2)<-c("Patient", "Count")
counts2
counts2$Patient<-str_extract(counts2$Patient, "[0-9]{4}")
counts2
counts2$log2<-log2(counts2$Count+1)
sort(counts2$log2, decreasing = FALSE)  # min 0 and max 11152 (20% max = >8920  // 20% min <2230)

upper20<-counts2%>%slice_max(Count, n=185)
upper20
lower20<-counts2%>%slice_min(Count, n=185)
lower20



lowerandupper20<-rbind(lower20, upper20)
view(lowerandupper20)
nrow(lowerandupper20)


### Arrange survival data set ###

Sur_df<-read_tsv("MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL.tsv")
Sur_df<-data.frame(Sur_df)
#Sur_df<-data.frame(MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL)
ncol(Sur_df)
Sur_df$Patient<-str_extract(Sur_df$PUBLIC_ID, "[0-9]{4}")
nrow(Sur_df)
Sur_df
Sur_df1<-Sur_df[1:1143, 1:3]
Sur_df1$Patient<-Sur_df$Patient
Sur_df1

topbottomsurv<-merge(lowerandupper20, Sur_df1, by ="Patient", all = FALSE)
nrow(topbottomsurv)
view(topbottomsurv)
sort(topbottomsurv$log2, decreasing = FALSE)

topbottomsurv<-topbottomsurv%>%mutate("Percentage" = case_when(log2<6.108524 ~ "Bottom", 
                                                               log2>9.294621 ~ "Top"))

topbottomsurv
nrow(topbottomsurv)
nrow(topbottomsurv%>%filter(Percentage == "Top"))
nrow(topbottomsurv%>%filter(Percentage == "Bottom"))
topbottomsurv<-topbottomsurv%>%mutate("Status" = case_when(deathdy>1 ~ 1)%>%replace_na(0)) 
## Death 1 - Alive 0

view(topbottomsurv)
topbottomsurv$Months_Alive<-topbottomsurv$lstalive/30
topbottomsurv


km <- with(topbottomsurv, Surv(Months_Alive, Status))
head(km,80)

km_fit <- survfit(Surv(Months_Alive, Status) ~ 1, data=topbottomsurv)
summary(km_fit)

ggsurvplot(km_fit, topbottomsurv,
           conf.int = TRUE,
           risk.table = TRUE,
           break.time.by=12,
           title="Overall Survival Rate",
           risk.table.height = 0.25,
           ggtheme = theme_bw())+
  xlab("Time (Days)")


by_percentage <- survfit(Surv(Months_Alive, Status) ~ Percentage, data=topbottomsurv)
summary(by_percentage)
by_percentage

ggsurvplot(by_percentage, topbottomsurv,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE,
           break.time.by=12,
           title="Overall Survival by CCR10 expression",
           risk.table.col = "strata",
           legend.labs =  c("Bottom 20%", "Top 20%"),   
           risk.table.height = 0.25, 
           ggtheme = theme_bw() )+
  xlab("Time (months)")


## Overall survival based on 50% top and 50% bottom

library(tidyverse)
library(tibble)
library(survminer)
library(survival)

################# Survival curve based on CCR10 expression (ENSG00000184451) ###################



## Read HtSeq data Survival based on 20% top and 20% bottom
raw<-read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
raw<-data.frame(raw)
#raw <- data.frame(Expression.Estimates...Gene.Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts)
ncol(raw)
raw<-raw%>%filter(Gene == "ENSG00000184451")
raw
counts <- as.matrix(raw[-1])
counts1<-data.frame(t(counts))
counts1
ncol(counts)



counts2<- tibble::rownames_to_column(counts1, "Patient")
colnames(counts2)<-c("Patient", "Count")
counts2
counts2$Patient<-str_extract(counts2$Patient, "[0-9]{4}")
counts2
counts2$log2<-log2(counts2$Count+1)
sort(counts2$log2, decreasing = FALSE)  # min 0 and max 11152 (20% max = >8920  // 20% min <2230)

upper50<-counts2%>%slice_max(Count, n=464)
view(upper50)
lower50<-counts2%>%slice_min(Count, n=464)
view(lower50)



lowerandupper50<-rbind(lower50, upper50)
view(lowerandupper50)
nrow(lowerandupper50)


### Arrange survival data set ###
Sur_df<-read_tsv("MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL.tsv")
Sur_df<-data.frame(Sur_df)

#Sur_df<-data.frame(MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL)
ncol(Sur_df)
Sur_df$Patient<-str_extract(Sur_df$PUBLIC_ID, "[0-9]{4}")
nrow(Sur_df)
Sur_df
Sur_df1<-Sur_df[1:1143, 1:3]
Sur_df1$Patient<-Sur_df$Patient
Sur_df1

topbottomsurv<-merge(lowerandupper50, Sur_df1, by ="Patient", all = FALSE)
nrow(topbottomsurv)
view(topbottomsurv)
sort(topbottomsurv$log2, decreasing = FALSE)

topbottomsurv<-topbottomsurv%>%mutate("Percentage" = case_when(log2<7.73 ~ "Bottom", 
                                                               log2>7.741467 ~ "Top"))

topbottomsurv
nrow(topbottomsurv)
nrow(topbottomsurv%>%filter(Percentage == "Top"))
nrow(topbottomsurv%>%filter(Percentage == "Bottom"))
topbottomsurv<-topbottomsurv%>%mutate("Status" = case_when(deathdy>1 ~ 1)%>%replace_na(0)) 
## Death 1 - Alive 0

view(topbottomsurv)
topbottomsurv$Months_Alive<-topbottomsurv$lstalive/30
topbottomsurv


km <- with(topbottomsurv, Surv(Months_Alive, Status))
head(km,80)

km_fit <- survfit(Surv(Months_Alive, Status) ~ 1, data=topbottomsurv)
summary(km_fit)

ggsurvplot(km_fit, topbottomsurv,
           conf.int = TRUE,
           risk.table = TRUE,
           break.time.by=12,
           title="Overall Survival Rate",
           risk.table.height = 0.25,
           ggtheme = theme_bw())+
  xlab("Time (Days)")


by_percentage <- survfit(Surv(Months_Alive, Status) ~ Percentage, data=topbottomsurv)
summary(by_percentage)
by_percentage

ggsurvplot(by_percentage, topbottomsurv,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE,
           break.time.by=12,
           title="Overall Survival by CCR10 expression",
           risk.table.col = "strata",
           legend.labs =  c("Bottom 50%", "Top 50%"),   
           risk.table.height = 0.25, 
           ggtheme = theme_bw() )+
  xlab("Time (months)")
