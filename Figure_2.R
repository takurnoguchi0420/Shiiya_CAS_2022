### Figure 2D: CD24 expression between EGFR-WT vs. EGFR-Mu###
# Source data: http://linkedomics.org/login.php

### Load a package###
library(tidyverse)

### Load the original data###
#SNV data
snv <- read.delim("Human__TCGA_LUAD__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt",sep = "\t")
#RNA-seq data
rna <- read.delim("Human__TCGA_LUAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", sep = "\t")
#Clinical data
clinical <- read.delim("Human__TCGA_LUAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", sep = "\t")

# Data formatting
snv_long <- 
snv %>% filter(attrib_name == "EGFR") %>% pivot_longer(cols = c(2:534)) %>% select(c(2,3))

rna_long <- 
rna %>% filter(attrib_name == "CD24") %>% pivot_longer(cols = c(2:516)) %>% select(c(2,3))

clinical_stage34 <- 
  clinical %>% 
  filter(attrib_name == "pathologic_stage") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(V1,"stageiii|stageiv")) 
  
colnames(clinical_stage34) <- c("name", "stage")

# Merge clinical,SNV, and RNA-seq data.
summary <- 
  left_join(clinical_stage34, snv_long, by = "name") %>% 
  left_join(rna_long, by = "name")
  
# Make the final format.
final <- 
summary %>% drop_na(value.y) %>%  
  drop_na(value.x) %>%  
  mutate(EGFR = 
    case_when(value.x == 0 ~ "neg",
                   TRUE ~ "pos")) %>% 
  select(c(1,4,5)) %>% 
  rename(cd24_log2_normalized=value.y) 
# Save the data in a csv format.
write.csv(final,"CD24_log2_EGFR_pos_vs_neg_TCGA_LUAD_stageiii_stageiv.csv",row.names = F)

### Visualization and stats ###
final %>% 
  ggplot(aes(x=EGFR, y=cd24_log2_normalized))+
  geom_violin()+
  geom_point(position = position_jitter(seed = 1, width = 0.2))+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 16, margin = margin(t = 8, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 16),
        axis.ticks.length.x = unit(0.3, "cm"),
        axis.ticks.length.y = unit(0.2, "cm"))

pos <- final %>% filter(EGFR=="pos")
pos <- pos$cd24_log2_normalized
neg <- final %>% filter(EGFR=="neg")
neg <- neg$cd24_log2_normalized

t.test(pos, neg)
