#setwd("~/DI-manuscript-rev/Pre-intervention/using-orig-2019-files/01_miRNA")
setwd("C:/Users/tanpf/OneDrive - A STAR/Sperm_diet_intervention_sncRNA/Manuscript_prep_sept2023/Figure2/")

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

counts <- read.delim("./data/MV1-piRNA-CPM-TMM.txt", header = T) # 17 samples
colnames(counts)
length(counts)
counts2 <- counts[,-1] # ensure to remove the extra column that come after importing the file
length(counts2)#2291 = ID+2290 piRNAs

rownames (counts2)
rownames (metrics)
identical(rownames(metrics), rownames(counts2))   # same order
################################
## Combine metrics and counts ## 
################################
data <- full_join(metrics, counts2, by = "ID")
length(data)#2342 (52+2290 (ID same))

###############Corr with Age###################################
outputa <- NULL
for(i in 53:length(data))
{
  rescor<- (cor.test(as.numeric(data$Age), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE))
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  pir <- colnames(data[i])
  new <- data.frame(pir, est, pval)
  outputa <- rbind(outputa, new)
 }
  
write.table(outputa, file="pirs-age-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputa$pval)
outputa <- outputa[i,]
outputa %>% head()

outputa[1:2,]

#            pir        est        pval
#941  piR_009572 -0.7156863 0.001235509
#1205 piR_013247  0.7034314 0.001629430

a <- ggplot(data, aes(x = Age, y = log(piR_009572,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(i)",x="Age", y= "piR_009572") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

b <- ggplot(data, aes(x = Age, y =  log(piR_013247,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(ii)",x="Age", y= "piR_013247") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b, ncol = 4, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_piRNA.png", width = 2200,height = 1200) 


###########Corr with BMI###########################################################################################################
outputb <- NULL
for(i in 53:length(data))
{
  
  rescor <- cor.test(as.numeric(data$BMI), log(as.numeric(data[,i],10)), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  pir <- colnames(data[i])
  new <- data.frame(pir, est, pval)
  outputb <- rbind(outputb, new)
}

write.table(outputb, file="pirs-bmi-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputb$pval)
outputb <- outputb[i,]
outputb %>% head()

outputb[1:2,]

#            pir        est         pval
#1056 piR_011527  0.7411049 0.0006642129
#1599 piR_017211 -0.6625772 0.0037511516

c <- ggplot(data, aes(x = BMI, y = log(piR_011527,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(iii)",x="BMI", y= "piR_011527") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

d <- ggplot(data, aes(x = BMI, y =  log(piR_017211,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(iv)",x="BMI", y= "piR_017211") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b,c,d, ncol = 4, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_piRNA.png", width = 2200,height = 1200) 


###########Corr with Sperm concentration###########################################################################################

outputc <- NULL

for(i in 53:length(data))
{
  
  rescor <- cor.test(as.numeric(data$Sperm.concentration..million.per.ml.), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  pir <- colnames(data[i])
  new <- data.frame(pir, est, pval)
  outputc <- rbind(outputc, new)
}

write.table(outputc, file="pirs-spconc-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


i <- order(outputc$pval)
outputc <- outputc[i,]
outputc %>% head()

outputc[1:2,]

#            pir        est         pval
#257  piR_002703 -0.8865728 2.138709e-06
#1844 piR_018904 -0.8681792 6.245260e-06

e <- ggplot(data, aes(x = Sperm.concentration..million.per.ml., y = log(piR_002703,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(v)",x="Sperm conc.", y= "piR_002703") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

f <- ggplot(data, aes(x = Sperm.concentration..million.per.ml., y =  log(piR_018904,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(vi)",x="Sperm conc.", y= "piR_018904") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b,c,d,e,f, ncol = 4, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_piRNA.png", width = 2200,height = 1200) 

###########Corr with Sperm motility#########################################
outputm <- NULL

for(i in 53:length(data))
{
  rescor <- cor.test(as.numeric(data$Percentage.of.sperm.motile....), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  pir <- colnames(data[i])
  new <- data.frame(pir, est, pval)
  outputm <- rbind(outputm, new)  
}

write.table(outputm, file="pirs-spmot-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputm$pval)
outputm <- outputm[i,]
outputm %>% head()

outputm[1:2,]

#            pir        est         pval
#463  piR_004895 -0.7506165 0.0005170996
#1334 piR_014553  0.6228520 0.0075689629

g <- ggplot(data, aes(x = Percentage.of.sperm.motile...., y = log(piR_004895,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(vii)",x="Sperm motility", y= "piR_004895") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

h <- ggplot(data, aes(x = Percentage.of.sperm.motile...., y =  log(piR_014553,10))) + 
     geom_point(colour = "salmon", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(viii)",x="Sperm motility", y= "piR_014553") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b,c,d,e,f,g,h, ncol = 4, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_piRNA.png", width = 2700,height = 1200) 
