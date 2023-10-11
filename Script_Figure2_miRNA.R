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

counts <- read.delim("./data/MV1-miRNA-CPM-TMM.txt", header = T) # 17 samples
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
length(data)

###############Corr with Age###################################

outputa <- NULL

for(i in 53:length(data))
{
  rescor<- (cor.test(as.numeric(data$Age), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE))
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  mir <- colnames(data[i])
  new <- data.frame(mir, est, pval)
  outputa <- rbind(outputa, new)
}

write.table(outputa, file="./output/mirs-age-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputa$pval)
outputa <- outputa[i,]
outputa %>% head()


outputa[1:2,]

#            mir        est       pval
#148  miR.23b.3p -0.6225490 0.00760691
#241 miR.3925.3p  0.6053922 0.01001777

a <- ggplot(data, aes(x = Age, y = log(miR.23b.3p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(i)",x="Age", y= "miR.23b.3p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

b <- ggplot(data, aes(x = Age, y =  log(miR.3925.3p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(ii)",x="Age", y= "miR.3925.3p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

ggarrange(a, b, ncol = 2, nrow = 1) %>% 
  ggexport(filename = "./Output/Figure2_miRNA.png", width = 2000,height = 850)

###########Corr with BMI########################################
outputb <- NULL
for(i in 53:length(data))
{
  rescor <- cor.test(as.numeric(data$BMI), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  mir <- colnames(data[i])
  new <- data.frame(mir, est, pval)
  outputb <- rbind(outputb, new)
}

write.table(outputb, file="mirs-bmi-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputb$pval)
outputb <- outputb[i,]
outputb %>% head()

outputb[1:2,]
#             mir        est        pval
#128   miR.204.5p -0.6588962 0.004019812
#155 miR.26a.2.3p  0.6184054 0.008141351

c <- ggplot(data, aes(x = BMI, y = log(miR.204.5p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(iii)",x="BMI", y= "miR.204.5p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

d <- ggplot(data, aes(x = BMI, y =  log(miR.26a.2.3p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(iv)",x="BMI", y= "miR.26a.2.3p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b,c,d, ncol = 2, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_miRNA.png", width = 2000,height = 850)



###########Corr with Sperm concentration#########################################
outputc <- NULL
for(i in 53:length(data))
{
  rescor <- cor.test(as.numeric(data$Sperm.concentration..million.per.ml.), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  mir <- colnames(data[i])
  new <- data.frame(mir, est, pval)
  outputc <- rbind(outputc, new)
}

write.table(outputc, file="mirs-spconc-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


i <- order(outputc$pval)
outputc <- outputc[i,]
outputc %>% head()

outputc[1:2,]
#          mir        est         pval
#180 miR.31.5p  0.8890253 1.828742e-06
#3   let.7b.5p -0.8669530 6.669552e-06

e <- ggplot(data, aes(x = Sperm.concentration..million.per.ml., y = log(miR.31.5p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(v)",x="Sperm conc.", y= "miR.31.5p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

f <- ggplot(data, aes(x = Sperm.concentration..million.per.ml., y =  log(let.7b.5p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(vi)",x="Sperm conc.", y= "let.7b.5p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b,c,d,e,f, ncol = 4, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_miRNA.png", width = 2000,height = 850)




###########Corr with Sperm motility#########################################
outputm <- NULL
for(i in 53:length(data))
{
  rescor <- cor.test(as.numeric(data$Percentage.of.sperm.motile....), log(as.numeric(data[,i]),10), method = "spearman", exact = FALSE)
  est <-as.numeric(rescor$estimate)
  pval <-as.numeric(rescor$p.value)
  mir <- colnames(data[i])
  new <- data.frame(mir, est, pval)
  outputm <- rbind(outputm, new)
}

write.table(outputm, file="mirs-spmot-corr1811.txt", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")

i <- order(outputm$pval)
outputm <- outputm[i,]
outputm %>% head()

outputm[1:2,]

#            mir        est       pval
#254    miR.4443  0.5761688 0.01549012
#127  miR.204.3p  0.5724833 0.01632072


g <- ggplot(data, aes(x = Percentage.of.sperm.motile...., y = log(miR.4443,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(vii)",x="Sperm motility ", y= "miR.4443") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))

h <- ggplot(data, aes(x = Percentage.of.sperm.motile...., y =  log(miR.204.3p,10))) + 
     geom_point(colour = "olivedrab3", size = 5) + 
     geom_smooth(method="lm", formula = y~x, colour="black") + 
     labs(title = "(viii)",x="Sperm motility", y= "miR.204.3p") + 
     stat_cor(method = "spearman", label.x.npc= "left", label.y.npc = "top", size= 12)   +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text() ,text = element_text(size = 50,color="black"), axis.text =element_text(color="black", size = 50))


ggarrange(a,b,c,d,e,f,g,h, ncol = 4, nrow = 2) %>% 
  ggexport(filename = "./Output/Figure2_miRNA.png", width = 2700,height = 1200) 
