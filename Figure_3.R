### Load packages ###
library(ggsci)
library(tidyverse)
library(readxl)

### Figure 3A: Fraction of cell types in the tumor microenvironment in NSCLC ###
# TCGA-LUDA clinical data source: http://linkedomics.org/login.php
# TCGA-LUDA deconvolution data source: https://github.com/BostonGene/Kassandra

# Load the clinical data.
clinical <- read.delim("/path/to/Human__TCGA_LUAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", sep = "\t")

# Select stage III and IV patients.
clinical_stage34 <- 
  clinical %>% 
  filter(attrib_name == "pathologic_stage") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(str_detect(V1,"stageiii|stageiv")) %>% 
  rename(patient_id = rowname, stage = V1) %>% 
  mutate(patient_id = str_replace_all(patient_id, "\\.", "-")) %>% 
  tibble()

# Load the Kassandra data.
kas <- read_tsv("/path/to/Kassandra/data/tcga-kassandra/LUAD_BG_deconvolution.tsv")

# Define a list of cell types to be analyzed.
cell_type_filtered <- c("B_cells",
                        "CD8_T_cells_PD1_high",
                        "CD8_T_cells_PD1_low",
                        "T_helpers",
                        "Tregs",
                        "NK_cells",
                        "Macrophages_M1",
                        "Macrophages_M2",
                        "Monocytes",
                        "Neutrophils",
                        "Endothelium",
                        "Fibroblasts",
                        "Other")

# Generate a tidy data for visualization.
# Cell fractions in StageIII/IV and recurrent patients will be analyzed.
# If you are unfamiliar with the TCGA barcoding method, please refer to:
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
kas_l <- 
  kas %>% 
  pivot_longer(cols = c(-1),values_to = "fraction") %>% 
  rename(cell_type = ...1) %>% 
  mutate(sample_type = str_extract(name, "(?<=^.{13,13})..")) %>% 
  mutate(suffix = str_extract(name, "(?<=^.{16,16}).*")) %>%  
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(patient_id = str_extract(name, "^.{12,12}")) %>% 
  filter(cell_type %in% cell_type_filtered) %>% 
  mutate(cell_type=fct_relevel(cell_type,cell_type_filtered)) %>% 
  filter(sample_type != "11") %>% #remove tissue normal
  filter(patient_id %in% clinical_stage34$patient_id|sample_type == "02")

# Change the data in a cross-table for PRISM.
kas_l_wide <- 
  kas_l %>% 
  select(c(1,3,6)) %>% 
  pivot_wider(names_from = cell_type,
              values_from = fraction) %>% 
  select(c(1,cell_type_filtered))

write_csv(kas_l_wide, file = "tcga-luad-StageiiiStageivRecurrent_kassandra.csv")

# Visualization

kas_l %>% 
  ggplot(aes(x=cell_type, y=fraction))+
  geom_boxplot(outlier.size = 1)+
  theme_bw()+
  ylab("Fraction")+
  xlab("Cell types")+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14,
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), 
                                   size = 12,
                                   angle = 20),
        axis.text.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0),
                                   size = 12),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        legend.text = element_text(size = 12))

### Figure 3B:Siglec-10 expression across immune cells###

#Zilionis, R. et al. Single-Cell Transcriptomics 
#of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations 
#across Individuals and Species. Immunity 50, 1317-1334.e10 (2019).
# Make first Table_S2_human_tumor_imm_minor.csv from the original Table_S2.xlsb.

# Load the original data.
dat <- read_csv("Table_S2_human_tumor_imm_minor.csv",
                col_names = TRUE)

# Siglec-10 expression across immune cell types
dat_siglec10 <- 
  dat %>% 
  filter(GENE == "SIGLEC10") %>% 
  pivot_longer(cols = -1, names_to = "subtypes") 

# Save the data in a csv format for PRISM.
write_csv(dat_siglec10,file = "siglec10_tils_scrnaseq.csv")

# Visualization
dat_siglec10 %>% 
  ggplot(aes(x=subtypes,y=value))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  xlab("Immune cell type")+
  ylab("Arbitrary expression value")+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14,
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), 
                                   size = 12,
                                   angle = 50),
        axis.text.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0),
                                   size = 12),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        legend.text = element_text(size = 12))
        
###Supplementary Figure 3: Sigle-10 and MRC1 expression in macrophage clusters###
dat %>% 
  filter(GENE %in% c("SIGLEC10","MRC1")) %>% 
  pivot_longer(cols = -1, names_to = "subtypes") %>% 
  filter(str_detect(subtypes,"hMø")) %>%
  mutate(value = if_else(X == "SIGLEC10", value*10, value)) %>% 
  rename(Gene = X) %>% 
  ggplot(aes(x=subtypes,y=value, fill = Gene))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  scale_fill_jco()+
  scale_y_continuous(
    
    # Features of the first axis
    name = "MRC1 expression",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*0.1, name="SIGLEC10 expression")
  )+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), angle = 30),
        axis.text.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        legend.text = element_text(size = 12))


ggsave("{file name}.pdf", 
       plot = last_plot(), 
       width = {number},
       height = {number}, 
       units = "in",
       dpi = 300)
