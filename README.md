# IAA in pneumonia

This is the code used for the main analyses in "Effect of the gut microbiota-derived tryptophan metabolite indole-3-acetic acid in pneumonia" (submitted). For questions: Bob Kullberg at r.f.j.kullberg@amsterdamumc.nl. 

A file containing mock clinical metadata is included to show the type of data, data structure and its definitions ('Mock data and definitions.xlsx'). Data protection regulations do not allow public sharing of individual patient data. Please refer to the data availability statement in the manuscript for information on how to gain data access. 

Analyses can be run in R/Rstudio (for installation, see: https://posit.co/download/rstudio-desktop/). Session information (including version numbers of all packages) are provided below. No earlier/later versions of packages and R/Rstudio were tested. On a normal desktop computer, run time for these analyses is < 10 minutes, although transcriptome analyses might take a few extra minutes. 


## Step 1 - Load libraries

```
library(tidyverse)
library(tableone)
library(rio)
library(readxl)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(survival)
library(survminer)
library(clusterProfiler)
```

## Step 2 - Tryptophan metabolites are altered during severe pneumonia and correlate with mortality

```
tryp <- tryp %>%
  mutate(mortalityd90 = if_else(mortalityd90 == "1", "death", "survivor", "control")) %>%
  mutate(mortalityd90 = fct_relevel(mortalityd90, "death", "survivor", "control"))
tryp.heatmap <- tryp[,c(4:26)] # Select metabolites
tryp.heatmap <- mutate_all(tryp.heatmap, function(x) as.numeric(as.character(x)))
tryp.heatmap <- mutate_all(tryp.heatmap, function(x) scale(x)) %>%
  dselect(`Tryptophane`, 
          `Indole-3-Acetic acid`,  `Indole-3-Aldehyde`, `Indole-3-Lactic acid`, `Indole-3-Sulfate`,  `sum_indole`, 
          `Serotonine`, `5-OH-Tryptophane`, `5-OH-Indole acetic acid`, `N-acetyl-serotonine`, `Melatonine`, `sum_5ht`, 
          `Kynurenine`, `Kynurenic acid`, `Picolinic acid`, `3-OH-Kynurenine`, `Quinolinic acid`, `Xanthurenic acid`, `3-OH-Anthranilic acid`, `sum_ido`, 
          `ratio_kyn_iaa`, `ratio_iaa_trp`, `ratio_kyn_trp`)
```
```
colours = colorRamp2(c(-1, 0, 2), c("#003366", "#EEEEEE", "#660033"))
row_ha = rowAnnotation(group = as.factor(tryp.all$group), mortality = as.factor(tryp.all$mortalityd90))

Heatmap(tryp.heatmap, 
        col = colours,
        row_title = "Patient", 
        cluster_columns = F,
        show_row_dend = F,
        cluster_rows = T, 
        left_annotation = row_ha,
        row_split = c(tryp$mortalityd90), 
        name = "Normalised value")
```

Reshape dataframe to enable Hedges' G calculation:
```
tryp.groups <- tryp.all %>%
  dselect(sample, group, 
          `Tryptophane`, 
          `Indole-3-Acetic acid`,  `Indole-3-Aldehyde`, `Indole-3-Lactic acid`, `Indole-3-Sulfate`,  `sum_indole`, 
          `Serotonine`, `5-OH-Tryptophane`, `5-OH-Indole acetic acid`, `N-acetyl-serotonine`, `Melatonine`, `sum_5ht`, 
          `Kynurenine`, `Kynurenic acid`, `Picolinic acid`, `3-OH-Kynurenine`, `Quinolinic acid`, `Xanthurenic acid`, `3-OH-Anthranilic acid`, `sum_ido`, 
          `ratio_kyn_iaa`, `ratio_iaa_trp`, `ratio_kyn_trp`) %>%
  reshape2::melt(id.vars = c("sample" ,"group"), 
                 variable.name = "metabolite")

effect_size_groups <- tryp.groups  %>%
  group_by(metabolite) %>%
  rstatix::cohens_d(value ~ group, hedges.correction = T, var.equal = F) %>% 
  as.data.frame()
```

Make the figure:
```
heatmap <- effect_size_groups %>% select(metabolite, effsize) %>%
  column_to_rownames("metabolite")
col_fun = colorRamp2(breaks=c(-1, 0, 1), c( "#660033", "#EEEEEE", "#003366"))

Heatmap(heatmap, 
        column_title = "Patients vs. Controls", 
        col = col_fun,
        show_column_names = F, 
        show_row_names = T, 
        row_names_side = "left",
        cluster_rows = F,
        cluster_columns = F,
        border=T,
        cluster_row_slices = F,
        heatmap_legend_param = list(at = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8),
                                    labels = c("<-0.8", "<-0.5", "<-0.2", ">-0.2 & <0.2", ">0.2", ">0.5", ">0.8"),
                                    title = "Hedges' g", color_bar = "discrete", fontsize = 8),
        rect_gp = gpar(col = "black", lwd = 0.5),
        border_gp = gpar(col = "black", lwd = 1))
```
Similar code was used to calculate Hedges' G for the comparison between survivors and non-survivors. 

Next, we assessed the association between IAA and 90-day mortality:
```
df <- df %>%
  mutate(mortalityd30 = if_else(mortalityd30 == "1", 1, 0)) %>%
  mutate(mortalityd90 = if_else(mortalityd90 == "1", 1, 0)) %>%
  mutate(mortality1y = if_else(mortality1y == "1", 1, 0))

summary(coxph(Surv(time_to_death90, mortalityd90) ~ log2(`Indole-3-Acetic acid`), data=df))
```

For visualisation and interpretation purposes, we stratified patients into three (equally large) tertiles based on their IAA level at ICU admission. 
```
df <- df %>%
  mutate(tertiles = ntile(as.numeric(`Indole-3-Acetic acid`), 3)) %>%
  mutate(tertiles = as.factor(if_else(tertiles == 1, 'Low IAA', if_else(tertiles == 2, 'Intermediate IAA', 'High IAA')))) %>%
  mutate(tertiles = fct_relevel(tertiles, 'Low IAA', 'Intermediate IAA', 'High IAA'))

summary(coxph(Surv(time_to_death90, mortalityd90) ~ tertiles, data=df)) 
```

Importantly, the association between IAA and mortality in our cohort of severe pneumonia patients was observed even following adjustment for sex, age, body mass index, disease severity (quantified using the Sequential Organ Failure Assessment [SOFA] score), causative pathogen and comorbidities (diabetes, malignancy, immunocompromised state, cardiovascular, renal, and respiratory disease) in a multivariable model. 
```
summary(coxph(Surv(time_to_death90, mortalityd90) ~ log2(`Indole-3-Acetic acid`) +
                Patient_gender + as.numeric(Patient_age_at_ICU_Admission)+ 
                bmi + as.numeric(SOFA_total) + 
                causative_pathogen +
                diabetes + malignancy + immunocompromised + cardiovascular+ renal+ respiratory,
              data=df))
```

Similar associations were found in sensitivity analysis with shorter (30 days) and longer (1 year) follow-up:
```
summary(coxph(Surv(time_to_death30, mortalityd30) ~ log2(`Indole-3-Acetic acid`), data=df)) 
summary(coxph(Surv(time_to_death1y, mortality1y) ~ log2(`Indole-3-Acetic acid`), data=df))
```

We visualised differences in IAA between survivors and non-survivors using a box / dot plot and survival curves:
```
df %>%
  mutate(mortalityd90 = as.factor(mortalityd90)) %>%
  ggplot(aes(x=mortalityd90, y=`Indole-3-Acetic acid`,  fill=mortalityd90)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.5) +
  geom_jitter(color = "black", pch = 21, alpha =.9, size = 3, width = 0.1)+
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("") +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  ylab("Plasma indole-3-acetic acid (ng/mL)") +
  scale_fill_manual(values=c("#e6ab02", "#d95f02")) + 
  stat_compare_means(size=4, label = "p.format", comparisons = list(c("0", "1")))

fit <- survfit(Surv(time_to_death90, mortalityd90) ~ tertiles, data=df)
iaa.surv <- ggsurvplot(fit, data = df,
                       axes.offset = F, 
                       censor = F, 
                       xlab = "Days since ICU admission", 
                       show.legend = F,
                       xlim = c(0,90),
                       legend = "bottom",
                       palette = c( "#44AA99", "#DDCC77", "#882255"))
```

Higher IAA was associated with fewer ventilator-free days:
```
cor.test(df$vfd90, log2(df$`Indole-3-Acetic acid`))
df %>%
  ggplot(aes(y = vfd90, x =  `Indole-3-Acetic acid`)) +
  geom_point(size=2.5, colour = "#000000") +
  geom_smooth(method="lm", se=F, fullrange=FALSE, span = 0.99,  level=0.95, colour = "black") +
  theme_bw() +
  scale_x_continuous(trans = "log10")
```

Higher IAA was associated with higer respiratory SOFA scores:
```
df %>%
  ggplot(aes(x=SOFA_Respiration_highest, y=`Indole-3-Acetic acid`,  fill=SOFA_Respiration_highest)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.5) +
  geom_jitter(color = "black", pch = 21, alpha =.9, size = 3, width = 0.1)+
  theme_bw() + 
  theme(legend.position = "none") +
  xlab("") +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  ylab("Plasma indole-3-acetic acid (ng/mL)") +
  xlab("SOFA - respiratory") +
  scale_fill_manual(values = c("#4197cc", "#DDCC77", "#660033")) +
  stat_compare_means(size=4, label = "p.format", comparisons = list(c("0-1", "2"), 
                                                                    c("2", "3-4"), 
                                                                    c("0-1", "3-4")))
```

IAA was associated with a lower probability of bacteremia:
```
df.bacterial <- df %>%
  filter(causative_pathogen_combined == "Bacterial")
summary(glm(as.numeric(positiveBC) ~ log2(`Indole-3-Acetic acid`), family=binomial, data=df.bacterial))

df.bacterial %>%
  ggplot(aes(x=positiveBC, y=`Indole-3-Acetic acid`,  fill=positiveBC)) + 
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.5) +
  geom_jitter(color = "black", pch = 21, alpha =.9, size = 3, width = 0.1)+
  theme_bw() + 
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  xlab("") +
  ylab("Plasma indole-3-acetic acid (ng/mL)") +
  scale_fill_manual(values=c("#4197cc", "#660033")) + 
  stat_compare_means(comparisons = list(c("0", "1")))
```


## Step 3 - Dynamics of tryptophan metabolites during experimental bacterial pneumonia

We used our well-established mouse model to assess the dynamics of tryptophan metabolites during severe bacterial pneumonia 
```
tryp.mice <- mutate_all(tryp.mice, function(x) as.numeric(as.character(x)))
metabolites <- colnames(tryp.mice)
res.lmer = data.frame(variable = metabolites,
                      coefficient=0, sd=0, pvalue=2, padj=2, stringsAsFactors = F)

for(i in range) {
  ct = summary(lmer(log10(tryp.mice[,i]) ~ tryp.mice$timepoint + (1|tryp.mice$mouse_id)))[["coefficients"]]
  pvalue = Anova(lmer(log10(tryp.mice[,i]) ~ tryp.mice$timepoint + (1|tryp.mice$mouse_id)))
  res.lmer[i, "variable"] = metabolites[i]
  res.lmer[i, "coefficient"] = ct[[2]]
  res.lmer[i, "sd"] = ct[[4]]
  res.lmer[i, "pvalue"] = pvalue$`Pr(>Chisq)`
} 
res.lmer$padj = p.adjust(res.lmer$pvalue, method = "BH")
```
```
res.lmer <- res.lmer %>%
  mutate(variable = fct_relevel(variable, "ratio_kyn_trp",  "ratio_iaa_trp", "ratio_kyn_iaa",
                                "sum_ido", "3-OH-Anthranilic acid", "Xanthurenic acid", "Quinolinic acid", "3-OH-Kynurenine",  "Picolinic acid", "Kynurenic acid","Kynurenine", 
                                "sum_5ht",  "Melatonine", "N-acetyl-serotonine", "5-OH-Indole acetic acid", "5-OH-Tryptophane", "Serotonine", 
                                "sum_indole", "Indole-3-Sulfate", "Indole-3-Lactic acid", "Indole-3-Aldehyde", "Indole-3-Acetic acid",
                                "Tryptophane")) %>%
  mutate(change = if_else(coefficient < 0, "Decrease", "Increase")) %>%
  mutate(change = if_else(pvalue < 0.05, change, "Not significant"))

res.lmer %>%
  ggplot(aes(x=coefficient, y=variable, colour = change)) +
  geom_errorbar(aes(xmin=coefficient-(2*sd), xmax=coefficient+(2*sd)), width=.1) +
  geom_line()+
  geom_point() +
  scale_x_continuous(limits = c(-0.175, 0.175)) +
  geom_vline(xintercept = 0.0, linetype = "dashed") +
  theme_bw() +
  xlab("Estimate ±2SD") +
  ylab("") +
  scale_colour_manual(values = c("#4096c3", "#882255", "#777777")) + 
  theme(legend.position = "bottom")
```
Similar code was used to assess the dynamics in bronchoalveolar lavage fluid

## Step 4 - IAA aggravates pulmonary damage and reduces dissemination during bacterial pneumonia
```
cfu <- read_excel("~/Mice_combined_data.xlsx", sheet = "CFU") %>%
  mutate(timegroup = paste(time, group, sep = "_")) %>%
  filter(time != "CH_24h") %>%
  mutate(organ = fct_relevel(organ, "Lung", "Blood", "Liver", "Spleen"))
```
```
colors <- c("#3194CD", "#8C2928") #Set colors
comparisons <- list(c("24h_control", "24h_IAA"), c("42h_control", "42h_IAA"))

cfu %>%
  filter(organ == "Lung") %>%
  ggplot(aes(x=timegroup, y=cfu,  fill=group)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.1)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + #mean_se
  stat_summary(fun.y=mean, geom="crossbar", color="black", width=0.2) +
  scale_fill_manual(values = colors) +
  theme_bw() + 
  scale_y_continuous(trans = scales::log10_trans(), breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", math_format(10^.x)), limits = c(1, NA)) + 
  theme(legend.position = "none") +
  facet_wrap(~., scales = "free_y", nrow=1) +
  ylab("CFU/ml") +
  xlab("") +
  stat_compare_means(comparisons = comparisons, size=4, label = "p.value")
```
Similar code was used to compare and visualise CFUs in blood, lung pathology scores and organ damage markers

## Step 5 - IAA impacts bacterial pneumonia through AhR
```
cfu <- read_excel("~/Mice_combined_data.xlsx", sheet = "CFU") %>%
  mutate(timegroup = paste(time, group, sep = "_")) %>%
  filter(time == "CH_24h") %>%
  mutate(timegroup = paste(time, group, sep = "_")) %>%
  mutate(organ = fct_relevel(organ, "Lung", "Blood", "Liver", "Spleen"))
```
```
colors <- c("#3194CD",  "#8C2928", "#4d9b9f") #Set colors
comparisons <- list(c("control", "IAA"), c("IAA", "IAA+CH-2213191"), 
                    c("control", "IAA+CH-2213191"))

cfu %>%
  mutate(group = fct_relevel(group, "control", "IAA",  "IAA+CH-2213191")) %>%
  filter(organ == "Lung") %>%
  ggplot(aes(x=group, y=cfu,  fill=group)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.1)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + #mean_se
  stat_summary(fun.y=mean, geom="crossbar", color="black", width=0.2) +
  scale_fill_manual(values = colors) +
  theme_bw() + 
  scale_y_continuous(trans = scales::log10_trans(), breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", math_format(10^.x)), limits = c(1, NA)) + 
  theme(legend.position = "none") +
  ylab("CFU/ml") +
  xlab("") +
  stat_compare_means(comparisons = comparisons, size=4, label = "p.value")
```
Similar code was used to compare and visualise CFUs in blood, lung pathology scores and organ damage markers


## Step 6 - IAA increases ROS release during in vitro K. pneumoniae infection and in patients
```
ros <- read_excel("~/ROS_production.xlsx") %>%
  reshape2::melt(id.vars = c("sample", "group", "stimulus", "combined"), variable.name = "time") %>%
  mutate(value = as.numeric(value)) %>%
  mutate(time = as.numeric(time))

ros.short <- ros %>%
  dplyr::group_by(group, stimulus, combined,  time) %>%
  summarise(mean = mean(value), 
            sem = (sd(value)/6^0.5))
```
```
colors <- c("#3194CD", "#8C2928") 
ros.short %>%
  ggplot(aes(x=((time-1)*2.5), y=mean, colour = group, fill = group, shape = stimulus)) +
  geom_smooth(span = .5, se=F)+
  scale_fill_manual(values = colors) +
  scale_colour_manual(values = colors) +
  theme_bw()+
  ylab("ROS (chemiluminesence)") +
  xlab("Time (minutes)") +
  geom_uperrorbar(aes(ymax = mean+sem)) +
  geom_point(shape = 21, colour = "black", size = 2) +
  theme(legend.position = "none")
```
```
auc <- ros %>%
  dplyr::select(sample, combined, time, value) %>%
  dplyr::group_by(sample, combined) %>%
  dplyr::summarise(auc = DescTools::AUC(time, value, 
                                        method = "trapezoid",
                                        absolutearea = T)) %>%
  dplyr::ungroup() %>%
  mutate(combined = fct_relevel(combined, "control_klebs", "IAA_klebs", "control_pma", "IAA_pma"))

comparisons <- list(c("control_klebs", "IAA_klebs"),
                    c("control_pma", "IAA_pma"))
colors <- c("#3194CD", "#8C2928", "#3194CD", "#8C2928") 

auc %>%
  ggplot(aes(x=combined, y = auc, fill = combined)) +
  geom_jitter(color = "black", pch = 21, alpha =.85, size = 3, width = 0.1)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) + #mean_se
  stat_summary(fun.y=mean, geom="crossbar", color="black", width=0.2) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_colour_manual(values = colors) +
  scale_y_continuous(trans = scales::log10_trans(), breaks = c(10^3, 10^4, 10^4.5, 10^5),
                     labels = scales::trans_format("log10", math_format(10^.x))) + 
  ylab("ROS (total chemiluminesence)") +
  xlab("") +
  stat_compare_means(comparisons = comparisons) +
  theme(legend.position = "none")
```

Finally, we assessed the relevance of the association between IAA and increased ROS in patients with severe pneumonia using whole blood transcriptomics data
```
library(reactome.db)
library(formulaic)
genenames <- colnames(transcriptome)[c(49479, 2:49387)]
genenames <- formulaic::add.backtick(genenames, include.backtick = "as.needed")
genes <- as.data.frame(transcriptome[,c(49479, 2:49387)])
genes <- mutate_all(genes, function(x) as.numeric(as.character(x)))
```
```
res.transcriptome = data.frame(variable = "IAA", 
                               genename=genenames, 
                               Estimate=0, Std=0, Tvalue=0, adjR=0, Pvalue=2, Padj=2, stringsAsFactors = F)

template <- "log2(`Indole-3-Acetic acid`) ~ "

for(i in 2:49387){
  f = genes[i]
  frm = as.formula(sprintf("%s + %s", template, f))
  model = lm(frm, genes)
  coef = data.frame(summary(model)$coefficients)
  colnames(coef) = c("Estimate", "Std", "Tvalue", "Pvalue")
  res.transcriptome[i, "genename"] = genenames[i]
  res.transcriptome[i, "Pvalue"] = coef$Pvalue[2]
  res.transcriptome[i, "Estimate"] = coef$Estimate[2]
  res.transcriptome[i, "Std"] = coef$Std[2]
  res.transcriptome[i, "Tvalue"] = coef$Tvalue[2]
  res.transcriptome[i, "R"] = summary(model)[["r.squared"]]
}

res.transcriptome <- res.transcriptome %>%
  filter(genename != "`Indole-3-Acetic acid`")
```
Reactome pathways related to ‘Immune System’ (R-HSA-168256) were assessed
```
xx <- as.list(reactomePATHID2EXTID)
if(length(xx) > 0){
  # Get the value of the first key
  xx[[1]]
  # Get the values for multiget for a few keys
  if(length(xx) >= 3){
    xx[1:3]
  }
}
entrez.immunesystem <- xx$`R-HSA-168256`

res <- res.transcriptome %>%
  left_join(genesymbols) %>%
  filter(`Entrez Gene` %in% entrez.immunesystem) %>%
  group_by(`Entrez Gene`) %>% 
  filter(Pvalue == min(Pvalue)) 

genelist <- res$Tvalue
names(genelist) <- res$`Entrez Gene` 
genelist <- na.omit(genelist) # omit any NA values 
genelist = sort(genelist, decreasing = TRUE) # sort the list in decreasing order (required for clusterProfiler)

library(ReactomePA)
gse <- gsePathway(geneList=genelist, 
                  maxGSSize = 1200,
                  pvalueCutoff = 0.05, 
                  verbose = T,
                  seed = T,
                  pAdjustMethod = "BH")
gsea <- as.data.frame(gse)  

gsea %>% 
  filter(p.adjust < 0.05) %>%
  mutate(Sign = if_else(NES > 0, "Positive", "Negative")) %>%
  ggplot() +
  geom_bar(aes(y= NES, x=reorder(Description, NES), fill = Sign), stat="identity", color= "black") +
  coord_flip()+
  scale_fill_manual(values = c("#44AA99", "#882255")) +
  ylab("Normalized Enrichment Score")+
  xlab("") +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = "none")
```
```
x2 <- enrichplot::pairwise_termsim(gse)
gseaplotros <- enrichplot::gseaplot2(x2, geneSetID = 141, title = x2$Description[141], subplots = 1:3)
gseaplotros <- ggplotify::as.ggplot(gseaplotros)
```


## Session info
```
R version 4.2.3 (2023-03-15)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.4

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ReactomePA_1.42.0     formulaic_0.0.8       reactome.db_1.82.0    AnnotationDbi_1.60.2 
 [5] IRanges_2.32.0        S4Vectors_0.36.2      Biobase_2.58.0        BiocGenerics_0.44.0  
 [9] clusterProfiler_4.6.2 survminer_0.4.9       ggpubr_0.6.0          survival_3.7-0       
[13] data.table_1.15.4     ComplexHeatmap_2.14.0 circlize_0.4.16       reshape2_1.4.4       
[17] readxl_1.4.3          rio_1.2.1             tableone_0.13.2       lubridate_1.9.3      
[21] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
[25] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.1        
[29] tidyverse_2.0.0       phyloseq_1.42.0       DOSE_3.24.2          

loaded via a namespace (and not attached):
  [1] utf8_1.2.4                  R.utils_2.12.3              tidyselect_1.2.1           
  [4] RSQLite_2.3.7               BiocParallel_1.32.6         Rtsne_0.17                 
  [7] scatterpie_0.2.3            munsell_0.5.1               codetools_0.2-20           
 [10] withr_3.0.1                 colorspace_2.1-1            GOSemSim_2.24.0            
 [13] knitr_1.48                  rstudioapi_0.16.0           DescTools_0.99.55          
 [16] ggsignif_0.6.4              MatrixGenerics_1.10.0       GenomeInfoDbData_1.2.9     
 [19] KMsurv_0.1-5                polyclip_1.10-7             bit64_4.0.5                
 [22] farver_2.1.2                rhdf5_2.42.1                downloader_0.4             
 [25] treeio_1.22.0               vctrs_0.6.5                 generics_0.1.3             
 [28] gson_0.1.0                  xfun_0.46                   timechange_0.3.0           
 [31] R6_2.5.1                    doParallel_1.0.17           GenomeInfoDb_1.34.9        
 [34] splitstackshape_1.4.8       clue_0.3-65                 graphlayouts_1.1.1         
 [37] locfit_1.5-9.10             gridGraphics_0.5-1          bitops_1.0-8               
 [40] rhdf5filters_1.10.1         microbiome_1.20.0           cachem_1.1.0               
 [43] fgsea_1.24.0                DelayedArray_0.24.0         scales_1.3.0               
 [46] ggraph_2.2.1                enrichplot_1.18.4           rootSolve_1.8.2.4          
 [49] gtable_0.3.5                lmom_3.0                    tidygraph_1.3.1            
 [52] rlang_1.1.4                 GlobalOptions_0.1.2         splines_4.2.3              
 [55] lazyeval_0.2.2              rstatix_0.7.2               broom_1.0.6                
 [58] abind_1.4-5                 backports_1.5.0             qvalue_2.30.0              
 [61] tools_4.2.3                 ggplotify_0.1.2             biomformat_1.26.0          
 [64] RColorBrewer_1.1-3          proxy_0.4-27                Rcpp_1.0.13                
 [67] plyr_1.8.9                  zlibbioc_1.44.0             RCurl_1.98-1.14            
 [70] GetoptLong_1.0.5            viridis_0.6.5               cowplot_1.1.3              
 [73] zoo_1.8-12                  SummarizedExperiment_1.28.0 ggrepel_0.9.5              
 [76] cluster_2.1.6               fs_1.6.4                    survey_4.4-2               
 [79] magrittr_2.0.3              mvtnorm_1.2-5               matrixStats_1.3.0          
 [82] patchwork_1.2.0             hms_1.1.3                   xtable_1.8-4               
 [85] HDO.db_0.99.1               XML_3.99-0.16.1             gridExtra_2.3              
 [88] shape_1.4.6.1               compiler_4.2.3              shadowtext_0.1.4           
 [91] crayon_1.5.3                R.oo_1.26.0                 ggfun_0.1.5                
 [94] mgcv_1.9-1                  tzdb_0.4.0                  geneplotter_1.76.0         
 [97] aplot_0.2.3                 expm_0.999-9                Exact_3.3                  
[100] DBI_1.2.3                   tweenr_2.0.3                rappdirs_0.3.3             
[103] MASS_7.3-60.0.1             boot_1.3-30                 Matrix_1.6-5               
[106] ade4_1.7-22                 car_3.1-2                   permute_0.9-7              
[109] cli_3.6.3                   mitools_2.4                 R.methodsS3_1.8.2          
[112] parallel_4.2.3              igraph_2.0.3                GenomicRanges_1.50.2       
[115] pkgconfig_2.0.3             km.ci_0.5-6                 foreach_1.5.2              
[118] ggtree_3.6.2                annotate_1.76.0             multtest_2.54.0            
[121] XVector_0.38.0              yulab.utils_0.1.5           digest_0.6.36              
[124] graph_1.76.0                vegan_2.6-6.1               Biostrings_2.66.0          
[127] cellranger_1.1.0            fastmatch_1.1-4             tidytree_0.4.6             
[130] survMisc_0.5.6              gld_2.6.6                   graphite_1.44.0            
[133] rjson_0.2.21                lifecycle_1.0.4             nlme_3.1-165               
[136] jsonlite_1.8.8              Rhdf5lib_1.20.0             carData_3.0-5              
[139] viridisLite_0.4.2           fansi_1.0.6                 pillar_1.9.0               
[142] lattice_0.22-6              KEGGREST_1.38.0             fastmap_1.2.0              
[145] httr_1.4.7                  GO.db_3.16.0                glue_1.7.0                 
[148] png_0.1-8                   iterators_1.0.14            bit_4.0.5                  
[151] ggforce_0.4.2               class_7.3-22                stringi_1.8.4              
[154] blob_1.2.4                  DESeq2_1.38.3               memoise_2.0.1              
[157] e1071_1.7-14                ape_5.8 
```
