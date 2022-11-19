### The source data preparation ###
# Download the source data: GSE57156_RAW.tar.
# Untar the tar file to extract *.CELL.gz files.
# Unzip the gz files to extract *.CELL files.

# Download the human surfanome dataset from http://wlab.ethz.ch/cspa
# and save it as human_surfanome.csv

### Load packages ###
library(tidyverse)
library(affy)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)  # equivalent to GPL570
library(ggrepel)
library(pheatmap)

# Normalize the signals using RMA
out_f <- "GSE57156_rma_expression.txt"
dat <- ReadAffy()
eset <- rma(dat) #"exprs(eset)" returns a expression matrix
write.exprs(eset, file = out_f) # write the RMA normalized expression file

# Load the normalized expression matrix
exp <- read.delim("GSE57156_rma_expression.txt", header = TRUE, sep = "\t")

# Data wrangling HCC2279 DMSOvsErlotinib
colnames(exp)
exp <- exp %>% select(c("X", 
                        "GSM1376328_36585_HCC2279_SPR999_HGU133P.CEL",
                        "GSM1376329_36586_HCC2279_SPR999_HGU133P.CEL",
                        "GSM1376330_36587_HCC2279_SPR999_HGU133P.CEL",
                        "GSM1376331_36588_HCC2279_SPR999_HGU133P.CEL"))
rownames(exp) <- exp[,1]
exp <- exp[,2:5]
new_colname <- c(paste("Ctrl_",1:2, sep = ""), paste("Treated_", 1:2, sep = ""))
colnames(exp) <- new_colname

## QC the data matrix ##
# Check data distributions
boxplot(exp, log="y", las=3)
# Check if a mean is set at 500 per sample after trimming off the lowest 2%
apply(exp, 2, median, na.rm=TRUE)

# Clustering & visualization
data.dist <- as.dist(1-cor(exp, method = "pearson")) # Pearson for log2-transformed values
out <- hclust(data.dist, method = "average") # "average" is preferred to "ward".
plot(out)

###########
# DGE
###########
in_f <- exp # a log2 normalized matrix file
out_f <- "GSE57156HCC2779_DGE_all.txt"
param_G1 <- 3
param_G2 <- 3
# Set up a design matrix
groups <- factor(c(rep("Ctrl",2), rep("Treated",2))) # Grouping
design <- model.matrix(~ 0 + groups)
colnames(design) <- c("Ctrl", "Treated")
contrast.matrix <- makeContrasts(TvsC = Treated - Ctrl, levels = design)
fit <- lmFit(in_f, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
out <- eBayes(fit2)

tT <- topTable(out, adjust.method = "BH", number = nrow(in_f)) # DGE
tT <- rownames_to_column(tT)
DGEtbl <- left_join(tT,rownames_to_column(in_f), by = c("rowname" = "rowname"))

#Annotation
annotation(gse) # returning GPL570
columns(hgu133plus2.db)
annot <- AnnotationDbi::select(hgu133plus2.db, 
                               keys = tT$rowname, # a vector of PROBRIDs to be matched
                               columns = c("SYMBOL", "ENTREZID"), # retrieve matching IDs 
                               keytype = "PROBEID") # define the keytype to be mathched ie, PROBEID in this case

DGEtbl <- left_join(DGEtbl, annot, by = c("rowname" = "PROBEID"))
DGEtbl <- DGEtbl[!duplicated(DGEtbl$rowname),]# Drop off duplicate rows in PROBEIDs

# Append surface protein information
sftbl <- read.delim("human_surfanome.csv", sep = ",") #http://wlab.ethz.ch/cspa
sftbl <- sftbl[,c("UP_Protein_name", "ENTREZ_gene_ID", "ENTREZ.gene.symbol" )] # trimming sftbl
sftbl[,2] <- sapply(sftbl[,2], as.character)#sftbl$ENTREZ_gene_ID: integer -> charater
DGEtbl <- left_join(DGEtbl, sftbl, by = c("ENTREZID" = "ENTREZ_gene_ID"))
DGEtbl <- DGEtbl[!duplicated(DGEtbl$rowname),]# Drop off duplicate rows in PROBEIDs
write.table(DGEtbl, file = "GSE57156_DGE_all.csv", sep = ",")

# Make a data.frame only for surface proteins
dat <- read.delim("GSE57156_DGE_all.csv", sep = ",")
cd_dat <- dplyr::filter(dat, str_detect(SYMBOL,"^CD\\d")) # rows with "CD??"
cspa_dat <- dat[!is.na(dat$UP_Protein_name),] # rows with surface protein descriptions
sf_dat <-  full_join(cspa_dat, cd_dat)

# Check NAs in sf_dat
sum(is.na(sf_dat$ENTREZ.gene.symbol))

# Reorder by adj.P.Val
sf_dat <- sf_dat %>% 
  dplyr::arrange(adj.P.Val)

# Save the summary tables
write.table(sf_dat, "GSE57156_DGE_sf.csv", sep = ",")

##########################################
## Data Analysis
##########################################
dat_all <- read.delim("GSE57156_DGE_all.csv", sep = ",")
dat_sf <- read.delim("GSE57156_DGE_sf.csv", sep = ",")

# Remove duplicated rows by SYMBOL 
dat_alls <- dat_all %>% tibble() %>% 
  distinct(SYMBOL, .keep_all = TRUE)

dat_sfs <- dat_sf %>% tibble() %>% 
  distinct(SYMBOL, .keep_all = TRUE) 

# Volcano plot for all genes
g1 <- ggplot(dat_alls,aes(logFC, -log10(adj.P.Val))) + 
  geom_point(size = 1, color = "grey50")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  xlim(-4.5, 4.5)+
  theme(legend.title=element_blank(),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Volcano plot for surface genes with labeling
dat_sfs$genelabels <- ""
dat_sfs$genelabels <- ifelse(dat_sfs$SYMBOL == "CD24"
                             ,TRUE,FALSE)

g2 <- ggplot(dat_sfs,aes(logFC, -log10(adj.P.Val))) + 
  geom_point(color = ifelse(dat_sfs$SYMBOL == "CD24", "red", "grey50"))+
  geom_text_repel(aes(logFC, -log10(adj.P.Val)),
                  label = ifelse(dat_sfs$genelabels == TRUE, 
                                 as.character(dat_sfs$SYMBOL),""), 
                  box.padding = unit(0.45, "lines"), hjust=2, size = 6) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  xlim(-3.2, 3.2)+
  theme_bw()+
  theme(legend.title=element_blank(),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# off labels and tick numbers
g2 <- ggplot(dat_sfs,aes(logFC, -log10(adj.P.Val))) + 
  geom_point(color = ifelse(dat_sfs$SYMBOL == "CD24", "red", "black"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  xlim(-3.2, 3.2)+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

ggsave(filename = "GSE57156_all_volcano.png", plot = g1, height=10, width=16, dpi = 300)
ggsave(filename = "GSE57156_sf_volcano_poster.png", plot = g2, height=10, width=16, dpi = 300)


# Heatmap
matrix <- dat_alls %>% tibble() %>% 
  select(8:11)
pheatmap(matrix,
         color = colorRampPalette(c("blue","white","red"))(100),
         border_color = "black",
         cellwidth = 20, 
         cellheight = NA,
         scale = "row", 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_method = "average",
         treeheight_row = 0,
         show_rownames = FALSE,
         fontsize = 12, 
         fontsize_row = 12, 
         fontsize_col = 12,
         filename = "GSE57156_heatmap_all.png") 
