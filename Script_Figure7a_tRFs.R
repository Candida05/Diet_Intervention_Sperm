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

metrics <- read.delim("./data/MV1-sample-factors-ed.txt", header=T, na.strings=c("Missing","N/A","NA","Not tested")) %>% 
  mutate(GroupVisit = factor( paste0(Group, ".", Visit.number) ),
         CoupleID  = factor( paste0("P", Couple.ID)) )

table(metrics$GroupVisit)
table(metrics$Visit.number)
table(metrics$CoupleID)

####################################
## Read in counts expression data ##
####################################

counts <- read.delim("./data/MV1-tRF-CPM-TMM-new.txt", header = T) # 17 samples
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

}
write.table(outputa, file="trfs-age-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputa$pval)
outputa <- outputa[i,]
outputa %>% head()

outputa[1:2,]

#                     trf        est         pval
#3     chr1.tir5.4.GlyCCC  0.7965686 0.0001297636
#111 chr6.trf5b.40.ValTAC -0.7573529 0.0004302161

a <- ggplot(data, aes(x = Age, y = log(chr1.tir5.4.GlyCCC,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(i)",x="Age (years)", y= "chr1.tir5.4.GlyCCC") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

b <- ggplot(data, aes(x = Age, y =  log(chr6.trf5b.40.ValTAC,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(ii)",x="Age (years)", y= "chr6.trf5b.40.ValTAC") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))


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
}

write.table(outputb, file="trfs-bmi-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputb$pval)
outputb <- outputb[i,]
outputb %>% head()

outputb[1:2,]

#                     trf        est       pval
#114 chr6.trf5b.76.LysTTT  0.5717796 0.01648322
#137  chr7.trf5c.7.CysGCA -0.5288348 0.02906308
c <- ggplot(data, aes(x = BMI, y = log(chr6.trf5b.76.LysTTT,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(iii)",x=expression("BMI(kg/m"^2*")"), y= "chr6.trf5b.76.LysTTT") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

d <- ggplot(data, aes(x = BMI, y =  log(chr7.trf5c.7.CysGCA,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(iv)",x=expression("BMI(kg/m"^2*")"), y= "chr7.trf5c.7.CysGCA") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))


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
}
write.table(outputc, file="trfs-spconc-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputc$pval)
outputc <- outputc[i,]
outputc %>% head()

outputc[1:2,]

#                       trf        est         pval
#2   chr1.tir5.137.Undet... -0.8068671 9.075647e-05
#138     chr9.tir5.7.HisGTG -0.7369713 7.381636e-04

e <- ggplot(data, aes(x = Sperm.concentration..million.per.ml., y = log(chr1.tir5.137.Undet...,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(v)",x="Sperm conc. (million/ml)", y= "chr1.tir5.137.Undet...") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

f <- ggplot(data, aes(x = Sperm.concentration..million.per.ml., y =  log(chr9.tir5.7.HisGTG,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(vi)",x="Sperm conc. (million/ml)", y= "chr9.tir5.7.HisGTG") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))


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
}
write.table(outputm, file="trfs-spmot-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputm$pval)
outputm <- outputm[i,]
outputm %>% head()

outputm[1:2,]

#                      trf        est       pval
#57  chr19.trf5b.13.ValCAC -0.5933679 0.01204328
#123  chr6.trf5c.37.ValAAC -0.5798543 0.01469315

g <- ggplot(data, aes(x = Percentage.of.sperm.motile...., y = log(chr19.trf5b.13.ValCAC,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(vii)",x="Sperm motility (%)", y= "chr19.trf5b.13.ValCAC") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

h <- ggplot(data, aes(x = Percentage.of.sperm.motile...., y =  log(chr6.trf5c.37.ValAAC,10))) + 
     geom_point(colour = "cyan") + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(viii)",x="Sperm motility (%)", y= "chr6.trf5c.37.ValAAC") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 5)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))


ggarrange(a,b,c,d,e,f,g,h, ncol = 4, nrow = 2) %>% 
ggexport(filename = "./Output/Figure2_tRFs.pdf", width = 12, height = 8)
