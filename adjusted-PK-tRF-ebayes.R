setwd("~/DI-manuscript-rev/DE-tRF-analysis")

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

counts <- read.delim("PK-tRF-expression-profile-orig.txt", row.names="TID", header = T) # Raw tRF expression profile for all 34 samples
colnames(counts)

table(rownames(metrics) %in% colnames(counts))   # all found

colnames(counts)
rownames (metrics)
identical(rownames(metrics), colnames(counts))   # same order

raw <- DGEList(counts=counts, samples=metrics) 
rm(counts, metrics)
dim(raw) #5000

############
## Filter ##
############

rs <- rowSums(cpm(raw) >= 1)
plot(table(rs))
keep <- which(rs >= 8)
raw <- raw[keep, , keep.lib.sizes=FALSE]  # 5000 ->  349
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
dupCor$consensus.correlation # 0.324

v <- voomWithQualityWeights(norm, design, block=norm$samples$CoupleID, correlation=dupCor$consensus, plot=T)#plot=T
dupCor <- duplicateCorrelation(v, design, block=norm$samples$CoupleID)
dupCor$consensus.correlation # 0.324

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
write.table(tt.adjusted, file="DE-PK-trfs.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t") #main result

tt <- topTable(efit, coef="intervention_over_control", n=Inf) %>% 
  rownames_to_column("trfs") %>% 
  separate("trfs", into=c("chr", "type", "subtype", "codon"), remove=FALSE) %>% 
  mutate(LOD=-log10(P.Value),
         sig=sign(logFC)*( abs(logFC) > log2(1.5) & LOD >= -log10(0.01) ), # higher confidence 0.01
         newlabel=paste(chr, type, subtype, codon, sep="-")) %>%  
  select(trfs, type, newlabel, codon, logFC, P.Value, LOD, sig)

library(ggrepel)

tiff("volcano-DE-PK-trfs.tiff", height = 20, width = 20, units = "cm", res=300) #main result

ggplot(tt, aes(x=logFC, y=LOD)) + 
  geom_point(data=subset(tt, sig==0), size=1) +
  geom_point(data=subset(tt, sig!=0), col="lightskyblue2", aes(shape=type), size=4) + 
  theme_classic(base_size = 18) + 
  xlab("log2 fold change") + ylab("-log10 p-value") + 
  xlim(-4, 4) + ylim(0,4)+
  geom_text_repel(data=subset(tt, sig!=0), aes(label=newlabel, fontface=2), #col=type removed
                  max.iter=10000, seed=123456, box.padding=0.5, nudge_y=0.3) +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), lty=2) +
  geom_hline(yintercept=-log10(0.01), lty=2)
dev.off()
