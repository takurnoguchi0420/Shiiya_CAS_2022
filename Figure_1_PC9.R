# Load a normalized matrix
library("GEOquery")
library(Biobase)
gse <- getGEO("GSE75308", 
           GSEMatrix = TRUE,
           getGPL = FALSE)
gse <- gse[[1]]
fvarLabels(gse) <- make.names(fvarLabels(gse))
x <- exprs(gse) #Get expression set from gse object
gse$data_processing #  "The data were normalised using quantile normalisation with IlluminaGUI in R"
write.table(x = x, file = "GSE75308.expression.matrix.csv",sep = ",") #export expression matrix in file (.csv format).

# Load the normalized matrix
y <- read.delim("GSE75308.expression.matrix.csv", header = TRUE, ",")
y <- y[,1:6] # select samples for the subsequent analysis

## QC the data matrix
# Check data distributions
boxplot(y, log="y", las=3)

# Check if a mean is set at 500 per sample after trimming off the lowest 2%
apply(y, 2, median, na.rm=TRUE)

# Check NAs
sum(is.na(y)) # returns "2"
y[is.na(y)] <- 1 #replace na with 1

# Check if the normalized matrix has been transformed by log2
qx <- as.numeric(quantile(y, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { y[which(y <= 0)] <- NaN
z <- log2(y) } # the original matrix is not normalized if the object z is generated. 

# Clustering and visualization
data.dist <- as.dist(1-cor(z, method = "pearson")) # Pearson for log2-transformed values
out <- hclust(data.dist, method = "average") # "average" is preferred to "ward".
plot(out) # visualization

# Check the GSE names corresponding to the conditions
sample <- rbind(gse$geo_accession, gse$title) # make a sample to condition table.
sample <- sample[,1:6] # select samples untreated and treated.

##############################
#####DGE
##############################
in_f <- z # a log2 normalized matrix file
out_f <- "GSE75308_DGE.txt"
param_G1 <- 3
param_G2 <- 3

# Set up a design matrix
colnames(z) <- c(paste("Ctrl_",1:param_G1, sep = ""), paste("Treated_", 1:param_G2, sep = ""))
groups <- factor(c(rep("Ctrl",param_G1), rep("Treated", param_G2))) # Grouping
design <- model.matrix(~ 0 + groups)
colnames(design) <- c("Ctrl", "Treated")
contrast.matrix <- makeContrasts(TvsC = Treated - Ctrl, levels = design)
library(limma)
fit <- lmFit(in_f, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
out <- eBayes(fit2)

tT <- topTable(out, adjust.method = "BH", number = nrow(in_f)) # DGE
tT <- rownames_to_column(tT)
DGEtbl <- left_join(tT,rownames_to_column(in_f), by = c("rowname" = "rowname"))

#Annotation
annotation(gse) # returning GPL10588
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("illuminaHumanv4.db") # equivalent to GPL10588
library(illuminaHumanv4.db)

annot <- AnnotationDbi::select(illuminaHumanv4.db, 
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
write.table(DGEtbl, file = "GSE75308_DGE_all.csv", sep = ",")

# Make a data.frame only for surface proteins
dat <- read.delim("GSE75308_DGE_all.csv", sep = ",")
cd_dat <- dplyr::filter(dat, str_detect(SYMBOL,"^CD\\d")) # rows with "CD??"
cspa_dat <- dat[!is.na(dat$UP_Protein_name),] # rows with surface protein descriptions
sf_dat <-  full_join(cspa_dat, cd_dat)
# Check NAs in sf_dat
sum(is.na(sf_dat$ENTREZ.gene.symbol))
# Reorder by adj.P.Val
sf_dat <- sf_dat %>% 
  dplyr::arrange(adj.P.Val)

# Save the summary tables
write.table(sf_dat, "GSE75308_DGE_sf.csv", sep = ",")

##########################################
## Data Analysis
##########################################
dat_all <- read.delim("GSE75308_DGE_all.csv", sep = ",")
dat_sf <- read.delim("GSE75308_DGE_sf.csv", sep = ",")
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
        
# off labels and tick numbers
g1 <- ggplot(dat_alls,aes(logFC, -log10(adj.P.Val))) + 
  geom_point(size = 1, color = "grey50")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  xlim(-4.5, 4.5)+
  theme(legend.title=element_blank(),
        axis.title  = element_blank(),
        axis.text = element_blank(),
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
  xlim(-4.5, 4.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        
# off labels and tick numbers
g2 <- ggplot(dat_sfs,aes(logFC, -log10(adj.P.Val))) + 
  geom_point(size = 2, color = ifelse(dat_sfs$SYMBOL == "CD24", "red", "black"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  xlim(-4.5, 4.5)+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1.2))

ggsave(filename = "GSE75308_all_volcano1.png", plot = g1, height=10, width=16, dpi = 300)
ggsave(filename = "GSE75308_sf_volcano1.png", plot = g2, height=10, width=8, dpi = 300)

#####
# heatmap 
#####
library(pheatmap)
matrix <- dat_alls %>% tibble() %>% 
  select(8:13)
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
         filename = "GSE75308_heatmap_all.png") 
