# Setting working directory #

setwd("C:/Users/sanad/Desktop/R-codes")

# Library used in the analysis #

library(limma)
library(edgeR)

# Reading the data into the RStudio #

rawdata1_1 <- read.delim("abundance1_1.tsv")
rawdata1_2 <- read.delim("abundance1_2.tsv")
rawdata1_3 <- read.delim("abundance1_3.tsv")
rawdata1_4 <- read.delim("abundance1_4.tsv")
rawdata1_5 <- read.delim("abundance1_5.tsv")
rawdata1_6 <- read.delim("abundance1_6.tsv")
rawdata2_1 <- read.delim("abundance2_1.tsv")
rawdata2_2 <- read.delim("abundance2_2.tsv")
rawdata2_3 <- read.delim("abundance2_3.tsv")
rawdata2_4 <- read.delim("abundance2_4.tsv")
rawdata2_5 <- read.delim("abundance2_5.tsv")
rawdata2_6 <- read.delim("abundance2_6.tsv")
rawdata3_1 <- read.delim("abundance3_1.tsv")
rawdata3_2 <- read.delim("abundance3_2.tsv")
rawdata3_3 <- read.delim("abundance3_3.tsv")
rawdata3_4 <- read.delim("abundance3_4.tsv")
rawdata3_5 <- read.delim("abundance3_5.tsv")
rawdata3_6 <- read.delim("abundance3_6.tsv")


# preparing the data for the analysis #

object_to_extract_rows <- rawdata1_1
count_for_all_conditions1 <- object_to_extract_rows[,-1] 
rownames(count_for_all_conditions1) <- object_to_extract_rows[,1] 


ctrl_L_1<-rawdata1_1$est_counts
ctrl_D_1<-rawdata1_2$est_counts
TRIG_L_1<-rawdata1_3$est_counts
TRIG_D_1<-rawdata1_4$est_counts
CLIMB_L_1<-rawdata1_5$est_counts
CLIMB_D_1<-rawdata1_6$est_counts
ctrl_L_2<-rawdata2_1$est_counts
ctrl_D_2<-rawdata2_2$est_counts
TRIG_L_2<-rawdata2_3$est_counts
TRIG_D_2<-rawdata2_4$est_counts
CLIMB_L_2<-rawdata2_5$est_counts
CLIMB_D_2<-rawdata2_6$est_counts
ctrl_L_3<-rawdata3_1$est_counts
ctrl_D_3<-rawdata3_2$est_counts
TRIG_L_3<-rawdata3_3$est_counts  
TRIG_D_3<-rawdata3_4$est_counts 
CLIMB_L_3<-rawdata3_5$est_counts 
CLIMB_D_3<-rawdata3_6$est_counts


count_for_all_conditions2 <- cbind(count_for_all_conditions1,ctrl_L_1, ctrl_L_2, ctrl_L_3, ctrl_D_1, ctrl_D_2, ctrl_D_3, TRIG_L_1, TRIG_L_2, TRIG_L_3, TRIG_D_1, TRIG_D_2, TRIG_D_3, CLIMB_L_1, CLIMB_L_2, CLIMB_L_3, CLIMB_D_1, CLIMB_D_2, CLIMB_D_3)
count_for_all_conditions <- count_for_all_conditions2[, -(1:4)]
dim(count_for_all_conditions)

# Normalisation and filtration

group <- c("ctrl_L","ctrl_L","ctrl_L","ctrl_D","ctrl_D","ctrl_D","TRIG_L","TRIG_L","TRIG_L","TRIG_D","TRIG_D","TRIG_D","CLIMB_L","CLIMB_L","CLIMB_L","CLIMB_D","CLIMB_D","CLIMB_D")
count_DGEList <-DGEList(counts = count_for_all_conditions, group = group) # creating the DGEList
keptGENES_normalized <- calcNormFactors(count_DGEList,method="TMM") # 
keptGENES_normalized <- estimateCommonDisp(keptGENES_normalized)
keptGENES_normalized <- estimateTagwiseDisp(keptGENES_normalized)
keptGENES_normalized.CPM <- cpm(keptGENES_normalized,log=T,normalized.lib.sizes=T)
keptGENES <-rowSums(count_DGEList$counts > 10 ) >2
counts.keep <- keptGENES_normalized.CPM[keptGENES,]
dim(counts.keep)

############################################################################
#quality control of the data # data exploration before batch elemination#

columns_color <- c("cyan3","brown3","blueviolet","cyan3","brown3","blueviolet", "cyan3","brown3","blueviolet", "cyan3","brown3","blueviolet", "cyan3","brown3","blueviolet", "cyan3","brown3","blueviolet")
# MDS plot #
plotMDS(counts.keep, main="MDS-plot of data before batch effect elemination",col=columns_color)
# cluster dendrogram according to the distance #
DEGLISTOBJECT_correlation = cor(counts.keep, method = "pearson")
dist <- as.dist(1- abs(DEGLISTOBJECT_correlation))
hc_tree <- hclust (dist, method="complete")
plot(hc_tree)


############################################################

#############################################################
# Differential expression analysis with linear model
# before batch elimination
#############################################################


counts.keep.logged = counts.keep  

type = c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
         "trig","trig","trig","trig","trig","trig",
         "climb","climb","climb","climb","climb","climb")
dl = c("light","light","light","dark","dark","dark",
       "light","light","light","dark","dark","dark",
       "light","light","light","dark","dark","dark")
rep = c("r1","r2","r3","r1","r2","r3","r1","r2","r3",
        "r1","r2","r3","r1","r2","r3","r1","r2","r3")
TS <- paste(type, dl,  sep=".") 
TS <- factor(TS, levels=c("ctrl.light","ctrl.dark","trig.light","trig.dark","climb.light","climb.dark"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
fit <- lmFit(counts.keep.logged, design)  


cont.matrix <- makeContrasts(
  LvsD_in_ctrl  =ctrl.light  - ctrl.dark,
  LvsD_in_trig  =trig.light  - trig.dark,
  LvsD_in_climb =climb.light - climb.dark,
  levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2) 

length(which(results[,"LvsD_in_ctrl"]!=0))
length(which(results[,"LvsD_in_trig"]!=0))
length(which(results[,"LvsD_in_climb"]!=0))


cont.matrix <- makeContrasts(
  ctrl_vs_trig_L  = ctrl.light - trig.light,
  ctrl_vs_trig_D  = ctrl.dark  - trig.dark,
  ctrl_vs_climb_L  = ctrl.light  - climb.light,
  ctrl_vs_climb_D  = ctrl.dark  - climb.dark,
  levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)

length(which(results[,"ctrl_vs_trig_L"]!=0))
length(which(results[,"ctrl_vs_trig_D"]!=0))
length(which(results[,"ctrl_vs_climb_L"]!=0))
length(which(results[,"ctrl_vs_climb_D"]!=0))


cont.matrix <- makeContrasts(
  
  trig_L_vs_climb_L  = trig.light  - climb.light,
  trig_D_vs_climb_D  = trig.dark  - climb.dark,
  levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)

length(which(results[,"trig_L_vs_climb_L"]!=0))
length(which(results[,"trig_D_vs_climb_D"]!=0))
