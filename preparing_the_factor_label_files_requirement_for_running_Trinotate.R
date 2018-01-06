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



#batch effect elemination and DE analysis

rep = c("r1","r2","r3","r1","r2","r3","r1","r2","r3",
        "r1","r2","r3","r1","r2","r3","r1","r2","r3")
data <- removeBatchEffect(counts.keep ,rep)

type = c("ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
         "trig","trig","trig","trig","trig","trig",
         "climb","climb","climb","climb","climb","climb")
dl = c("light","light","light","dark","dark","dark",
       "light","light","light","dark","dark","dark",
       "light","light","light","dark","dark","dark")

TS <- paste(type, dl, sep=".")
TS <- factor(TS, levels=c("ctrl.light","ctrl.dark","trig.light","trig.dark","climb.light","climb.dark"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
fit <- lmFit(data,design)

cont.matrix <- makeContrasts(
  ctrl_vs_trig_L  = ctrl.light - trig.light,
  ctrl_vs_trig_D  = ctrl.dark  - trig.dark,
  ctrl_vs_climb_L  = ctrl.light  - climb.light,
  ctrl_vs_climb_D  = ctrl.dark  - climb.dark,
  levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)


# checking the length of the list of DEG #
length(which(results[,"ctrl_vs_trig_L"]!=0))
length(which(results[,"ctrl_vs_trig_D"]!=0))
length(which(results[,"ctrl_vs_climb_L"]!=0))
length(which(results[,"ctrl_vs_climb_D"]!=0))


##########################################################################

# creating the factor label files # recuirement for running Trinotate program for transcriptome. annotation

##########################################################################

compString = "ctrl_vs_trig_L"
TXTFileName = paste(compString,".txt")
write.table(which(results[,"ctrl_vs_trig_L"]!=0), TXTFileName, col.names = TRUE, row.names = TRUE, sep = "\t")

myfile = read.table("ctrl_vs_trig_L .txt", header = TRUE)
#rownames(myfile)

factorlabel = rep("Effectofinfection_in_light_trig", 9743)  # annotation done with factorlabel_name "EffectofLight".
genesnames = rownames(myfile)
ctrl_vs_trig_L = data.frame(factorlabel,genesnames)
write.table(ctrl_vs_trig_L, file = "ctrl_vs_trig_L_F.txt", col.names = FALSE, row.names = FALSE, sep = "\t")

############################################################################

compString = "ctrl_vs_trig_D"
TXTFileName = paste(compString,".txt")
write.table(which(results[,"ctrl_vs_trig_D"]!=0), TXTFileName, col.names = TRUE, row.names = TRUE, sep = "\t")

myfile = read.table("ctrl_vs_trig_D .txt", header = TRUE)
#rownames(myfile)

factorlabel = rep("Effectofinfection_in_dark_trig", 9723)  # annotation done with factorlabel_name "EffectofLight".
genesnames = rownames(myfile)
ctrl_vs_trig_D = data.frame(factorlabel,genesnames)
write.table(ctrl_vs_trig_D, file = "ctrl_vs_trig_D_F.txt", col.names = FALSE, row.names = FALSE, sep = "\t")

##########################################################################


compString = "ctrl_vs_climb_L"
TXTFileName = paste(compString,".txt")
write.table(which(results[,"ctrl_vs_climb_L"]!=0), TXTFileName, col.names = TRUE, row.names = TRUE, sep = "\t")

myfile = read.table("ctrl_vs_climb_L .txt", header = TRUE)
#rownames(myfile)

factorlabel = rep("Effectofinfection_in_light_climb", 9028)  # annotation done with factorlabel_name "EffectofLight".
genesnames = rownames(myfile)
ctrl_vs_climb_L = data.frame(factorlabel,genesnames)
write.table(ctrl_vs_climb_L, file = "ctrl_vs_climb_L_F.txt", col.names = FALSE, row.names = FALSE, sep = "\t")

############################################################################

compString = "ctrl_vs_climb_D"
TXTFileName = paste(compString,".txt")
write.table(which(results[,"ctrl_vs_climb_D"]!=0), TXTFileName, col.names = TRUE, row.names = TRUE, sep = "\t")

myfile = read.table("ctrl_vs_climb_D .txt", header = TRUE)
#rownames(myfile)

factorlabel = rep("Effectofinfection_in_dark_climb", 9255)  # annotation done with factorlabel_name "EffectofLight".
genesnames = rownames(myfile)
ctrl_vs_climb_D = data.frame(factorlabel,genesnames)
write.table(ctrl_vs_climb_D, file = "ctrl_vs_climb_D_F.txt", col.names = FALSE, row.names = FALSE, sep = "\t")

###################################################################################################################################

