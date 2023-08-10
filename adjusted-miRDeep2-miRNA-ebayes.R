setwd("~/DI-manuscript-rev/DE-miRNA-analysis")

library(edgeR)
library(limma)
library(tidyverse)
library(reshape2)
library(gridExtra)

rm(list=ls())

###############
## Phenotype ##
###############

metrics <- read.delim("All-sample-factors-ed.txt", row.names=1, header=T, na.strings=c("Missing","N/A","NA","Not tested")) %>% #factors for all the samples
mutate(GroupVisit = factor( paste0(Group, ".", Visit.number) ),
CoupleID  = factor( paste0("P", Couple.ID)) )

table(metrics$GroupVisit)
table(metrics$Visit.number)
table(metrics$CoupleID)

##########################
## Read in and annotate ##
##########################

counts <- read.delim("miRDeep2-miRNA-expression-profile-orig.txt", row.names="MID", header = T) # Raw miRDeep2 miRNA expression profile for all 34 samples
colnames(counts)

table(rownames(metrics) %in% colnames(counts))   # all found

colnames(counts)
rownames (metrics)
identical(rownames(metrics), colnames(counts))   # same order


raw <- DGEList(counts=counts, samples=metrics) 
rm(counts, metrics)
dim(raw) #2654

############
## Filter ##
############

rs <- rowSums(cpm(raw) >= 1)
plot(table(rs))
keep <- which(rs >= 8)
raw <- raw[keep, , keep.lib.sizes=FALSE]  # 2654 -> 1179 
dim(raw)
rm(rs, keep)

###################
## Normalization ##
###################

norm <- calcNormFactors(raw, method="TMM")
rm(raw)
norm$samples
norm.cpm.vals = cpm(norm)
norm.lcpm.vals = cpm(norm, log=T)

#########
## DGE ##
#########

table(norm$samples$GroupVisit)

design <- model.matrix( ~ 0 + GroupVisit + Age + BMI + Sperm.concentration..million.per.ml. + Percentage.of.sperm.motile...., data=norm$samples)
colnames(design) <- gsub("GroupVisit", "", colnames(design))

v <- voomWithQualityWeights(norm, design, plot=T) #plot=T
dupCor <- duplicateCorrelation(v, design, block=norm$samples$CoupleID)
dupCor$consensus.correlation # 0.1939

v <- voomWithQualityWeights(norm, design, block=norm$samples$CoupleID, correlation=dupCor$consensus, plot=T)#plot=T
dupCor <- duplicateCorrelation(v, design, block=norm$samples$CoupleID)
dupCor$consensus.correlation # 0.1942

v.cpm.vals = cpm(v)

fit <- lmFit(v, design, block=norm$samples$CoupleID, correlation=dupCor$consensus)

cm <- makeContrasts(
  interventionEffect = (intervention.MV4 - intervention.MV1),
  controlEffect = (control.MV4 - control.MV1), 
  intervention_over_control = (intervention.MV4 - intervention.MV1) - (control.MV4 - control.MV1),
  levels=design)

fit2 <- contrasts.fit(fit, cm)

efit <- eBayes(fit2, trend=TRUE)

tt.adjusted <- topTable(efit, coef="intervention_over_control", n=Inf)
write.table(tt.adjusted, file="DE-miRDeep2-miRNAs.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t") #main result

tt <- topTable(efit, coef="intervention_over_control", n=Inf) %>% 
  rownames_to_column("mirs") %>% 
  mutate(LOD=-log10(P.Value), sig=sign(logFC)*( abs(logFC) > log2(1.5) & LOD >= -log10(0.01) )) %>% # higher confidence 0.01
  select(mirs, logFC, P.Value, LOD, sig)

library(ggrepel)

tiff("volcano-DE-miRDeep2-miRNAs.tiff", height = 20, width = 20, units = "cm", res=300) #main result

ggplot(tt, aes(x=logFC, y=LOD)) + 
  geom_point(data=subset(tt, sig==0), size=1) +
  geom_point(data=subset(tt, sig!=0), col="olivedrab3", size=2) + 
  theme_classic(base_size = 18) + 
  xlab("log2 fold change") + ylab("-log10 p-value") + 
  xlim(-10, 10) + ylim(0,5)+
  geom_text_repel(data=subset(tt, sig!=0), aes(label=mirs, fontface=2), 
                  max.iter=10000, seed=123456, box.padding=0.5, nudge_y=0.2) +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), lty=2) +
  geom_hline(yintercept=-log10(0.01), lty=2)
dev.off()


