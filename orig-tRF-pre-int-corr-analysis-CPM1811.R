setwd("~/DI-manuscript-rev/Pre-intervention/using-orig-2019-files/02_tRF")
library(edgeR)
library(limma)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(dplyr)
rm(list=ls())

###############
## Phenotype ##
###############

metrics <- read.delim("MV1-sample-factors-ed.txt", header=T, na.strings=c("Missing","N/A","NA","Not tested")) %>% 
  mutate(GroupVisit = factor( paste0(Group, ".", Visit.number) ),
         CoupleID  = factor( paste0("P", Couple.ID)) )

table(metrics$GroupVisit)
table(metrics$Visit.number)
table(metrics$CoupleID)

####################################
## Read in counts expression data ##
####################################

counts <- read.delim("MV1-tRF-CPM-TMM-new.txt", header = T) # 17 samples
colnames(counts)
length(counts)


table(rownames(metrics) %in% rownames(counts))   # all found

rownames (counts)
rownames (metrics)
identical(rownames(metrics), rownames(counts))   # same order
################################
## Combine metrics and counts ## 
################################
data <- full_join(metrics, counts, by = "ID")
length(data) #195
###############Corr with Age###################################
outputa <- NULL
for(i in 53:length(data))
{
  rescor<- (cor.test(as.numeric(data$Age), as.numeric(data[,i]), method = "spearman", exact = FALSE))
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  trf <- colnames(data[i])
  new <- data.frame(trf, est, pval)
  outputa <- rbind(outputa, new)

 if (rescor$p.value < 0.002)# get 6 at 0.002 cut-off 
  {
   print (paste(colnames(data)[i], " cor:", rescor$estimate, " p=value:", rescor$p.value))
   
    tiff(paste(colnames(data)[i],"Age.tiff"), height = 20, width = 20, units = "cm", res=300)
    p <- (ggplot(data, aes(x=Age,  y=data[,i])) + geom_point(colour="lightskyblue2", size =3 ) + geom_smooth(method="lm", colour="black")) + labs(x="Age", y=paste(colnames(data)[i],"(CPM)")) + theme_classic() + theme(axis.text = element_text(color = "black", size = 14), axis.title = element_text(color = "black", size = 16))
    print(p+(stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size=4)))   
    dev.off()
  }
}
write.table(outputa, file="trfs-age-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
###########Corr with BMI#######################################3##
outputb <- NULL
for(i in 53:length(data))
{
  
  rescor <- cor.test(as.numeric(data$BMI), as.numeric(data[,i]), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  trf <- colnames(data[i])
  new <- data.frame(trf, est, pval)
  outputb <- rbind(outputb, new)
  
  
  if (rescor$p.value < 0.03)# get 2 at 0.03 cut-off
  {
   print (paste(colnames(data)[i], " cor:", rescor$estimate, " p=value:", rescor$p.value))
   tiff(paste(colnames(data)[i],"BMI.tiff"), height = 20, width = 20, units = "cm", res=300)
   p <- (ggplot(data, aes(x=BMI, y=data[,i])) + geom_point(colour="lightskyblue2", size =3 ) + geom_smooth(method="lm", colour="black")) + labs(x="BMI", y=paste(colnames(data)[i],"(CPM)")) + theme_classic() + theme(axis.text = element_text(color = "black", size = 14), axis.title = element_text(color = "black", size = 16))
    print(p+(stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size=4)))   
    dev.off()
  }

}
write.table(outputb, file="trfs-bmi-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")
###########Corr with Sperm concentration#########################################
outputc <- NULL
for(i in 53:length(data))
{
  rescor <- cor.test(as.numeric(data$Sperm.concentration..million.per.ml.), as.numeric(data[,i]), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  trf <- colnames(data[i])
  new <- data.frame(trf, est, pval)
  outputc <- rbind(outputc, new)
  
  if (rescor$p.value < 0.002)#get 6 at 0.002 cut-off
  {
    print (paste(colnames(data)[i], " cor:", rescor$estimate, " p=value:", rescor$p.value))
    tiff(paste(colnames(data)[i],"Sp-conc.tiff"), height = 20, width = 20, units = "cm", res=300)
    p <- (ggplot(data, aes(x=Sperm.concentration..million.per.ml., y=data[,i])) + geom_point(colour="lightskyblue2", size =3 ) + geom_smooth(method="lm", colour="black")) + labs(x="Sperm concentration million per ml", y=paste(colnames(data)[i],"(CPM)")) + theme_classic() + theme(axis.text = element_text(color = "black", size = 14), axis.title = element_text(color = "black", size = 16)) 
    print(p+(stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size=4)))   
    dev.off()
  }
  
}
write.table(outputc, file="trfs-spconc-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

###########Corr with Sperm motility#########################################
outputm <- NULL
for(i in 53:length(data))
{
  rescor <- cor.test(as.numeric(data$Percentage.of.sperm.motile....), as.numeric(data[,i]), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  trf <- colnames(data[i])
  new <- data.frame(trf, est, pval)
  outputm <- rbind(outputm, new)
  
  if (rescor$p.value < 0.03)# get 3 at 0.03 cut-off
  {
    print (paste(colnames(data)[i], " cor:", rescor$estimate, " p=value:", rescor$p.value))
    tiff(paste(colnames(data)[i],"Sp-mot.tiff"), height = 20, width = 20, units = "cm", res=300)
    p <- (ggplot(data, aes(x=Percentage.of.sperm.motile...., y=data[,i])) + geom_point(colour="lightskyblue2", size =3 ) + geom_smooth(method="lm", colour="black")) + labs(x="Percentage of sperm motile", y=paste(colnames(data)[i],"(CPM)") ) + theme_classic() + theme(axis.text = element_text(color = "black", size = 14), axis.title = element_text(color = "black", size = 16))
    print(p+(stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size=4)))   
    dev.off()
  }
  
}
write.table(outputm, file="trfs-spmot-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

