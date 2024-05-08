rm(list = ls())

# LIBRARY
pacman::p_load(tidyverse,ggplot2, ggpubr, rstatix)

# SETPATH

dir <- "~/Sperm_diet_intervention_sncRNA/Manuscript_prep_sept2023/Figure3abc"

# READ DATA

sperm17s <- read.csv(paste0(dir,"/data/Data_SpermDietIntervention_17s.csv"), header = T)
head(sperm17s); dim(sperm17s)

sperm17s_mv1 <- subset(sperm17s, sperm17s$Visit == "MV1")
sperm17s_mv1 %>% head()
sperm17s_mv1 %>% names()

colnames(sperm17s_mv1)[6:51] <- paste0(colnames(sperm17s_mv1)[6:51], "_MV1")

sperm17s_mv4 <- subset(sperm17s, sperm17s$Visit == "MV4")
sperm17s_mv4 %>% head()
sperm17s_mv4 %>% names()

colnames(sperm17s_mv4)[6:51] <- paste0(colnames(sperm17s_mv4)[6:51], "_MV4")

# CHANGE IN VAR

table(sperm17s_mv1$CoupleID == sperm17s_mv4$CoupleID)

# CHANGE IN CONC. VITD IN SERUM

sperm17s <- sperm17s_mv1[,c(1:5,31)]
colnames(sperm17s)[6:6] <- "Group"
sperm17s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV1 <- sperm17s_mv1$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV1
sperm17s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV4 <- sperm17s_mv4$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV4 
sperm17s$Change_in_VitD_Conc_Serum_nmol_L <- sperm17s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV4 - sperm17s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV1

sperm17s %>% head()
sperm17s$EPA_in_RBC_percent_MV1 <- sperm17s_mv1$EPA_in_RBC_._MV1
sperm17s$EPA_in_RBC_percent_MV4 <- sperm17s_mv4$EPA_in_RBC_._MV4
sperm17s$Change_in_percent_EPA_in_RBC <- sperm17s$EPA_in_RBC_percent_MV4 - sperm17s$EPA_in_RBC_percent_MV1

sperm17s %>% head()
sperm17s$DHA_in_RBC_percent_MV1 <- sperm17s_mv1$DHA_in_RBC_._MV1
sperm17s$DHA_in_RBC_percent_MV4 <- sperm17s_mv4$DHA_in_RBC_._MV4
sperm17s$Change_in_percent_DHA_in_RBC <- sperm17s$DHA_in_RBC_percent_MV4 - sperm17s$DHA_in_RBC_percent_MV1

sperm17s %>% head()
sperm17s$EPA_in_Seminal_Fluid_MV1 <- sperm17s_mv1$EPA_in_Seminal_Fluid_._MV1
sperm17s$EPA_in_Seminal_Fluid_MV4 <- sperm17s_mv4$EPA_in_Seminal_Fluid_._MV4
sperm17s$Change_in_percent_EPA_in_Seminal_Fluid <- sperm17s$EPA_in_Seminal_Fluid_MV4 - sperm17s$EPA_in_Seminal_Fluid_MV1

sperm17s %>% head()
sperm17s$DHA_in_Seminal_Fluid_MV1 <- sperm17s_mv1$DHA_in_Seminal_Fluid_._MV1
sperm17s$DHA_in_Seminal_Fluid_MV4 <- sperm17s_mv4$DHA_in_Seminal_Fluid_._MV4
sperm17s$Change_in_percent_DHA_in_Seminal_Fluid <- sperm17s$DHA_in_Seminal_Fluid_MV4 - sperm17s$DHA_in_Seminal_Fluid_MV1

sperm17s %>% head()
sperm17s$EPA_in_Sperm_MV1 <- sperm17s_mv1$EPA_in_Sperm._MV1
sperm17s$EPA_in_Sperm_MV4 <- sperm17s_mv4$EPA_in_Sperm._MV4
sperm17s$Change_in_percent_EPA_in_Sperm <- sperm17s$EPA_in_Sperm_MV4 - sperm17s$EPA_in_Sperm_MV1

sperm17s %>% head()
sperm17s$DHA_in_Sperm_MV1 <- sperm17s_mv1$DHA_in_Sperm_._MV1
sperm17s$DHA_in_Sperm_MV4 <- sperm17s_mv4$DHA_in_Sperm_._MV4
sperm17s$Change_in_percent_DHA_in_Sperm <- sperm17s$DHA_in_Sperm_MV4 - sperm17s$DHA_in_Sperm_MV1


sperm17s$Group <- gsub("placebo", "Control", sperm17s$Group)
sperm17s$Group <- gsub("treatment", "Intervention", sperm17s$Group)

###

sperm102s <- read.csv(paste0(dir,"/data/Data_SpermDietIntervention_102s.csv"), header = T)
head(sperm102s); dim(sperm102s)

sperm102s_mv1 <- subset(sperm102s, sperm102s$Visit.number == "2")
sperm102s_mv1 %>% head()
sperm102s_mv1 %>% names()

colnames(sperm102s_mv1)[5:50] <- paste0(colnames(sperm102s_mv1)[5:50], "_MV1")

sperm102s_mv4 <- subset(sperm102s, sperm102s$Visit.number == "4")
sperm102s_mv4 %>% head()
sperm102s_mv4 %>% names()

colnames(sperm102s_mv4)[5:50] <- paste0(colnames(sperm102s_mv4)[5:50], "_MV4")


# PLOT BOXPLOT
table(sperm102s_mv1$Couple.ID == sperm102s_mv4$Couple.ID)

# CHANGE IN CONC. VITD IN SERUM

sperm102s <- sperm102s_mv1[,c(1:4,30)]
colnames(sperm102s)[5:5] <- "Group"
sperm102s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV1 <- sperm102s_mv1$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV1
sperm102s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV4 <- sperm102s_mv4$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV4 
sperm102s$Change_in_VitD_Conc_Serum_nmol_L <- sperm102s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV4 - sperm102s$Total_concentration_of_Vitamin_D_in_blood_serum_nmol_L_MV1

sperm102s %>% head()
sperm102s$EPA_in_RBC_percent_MV1 <- sperm102s_mv1$EPA_in_RBC_._MV1
sperm102s$EPA_in_RBC_percent_MV4 <- sperm102s_mv4$EPA_in_RBC_._MV4
sperm102s$Change_in_percent_EPA_in_RBC <- sperm102s$EPA_in_RBC_percent_MV4 - sperm102s$EPA_in_RBC_percent_MV1

sperm102s %>% head()
sperm102s$DHA_in_RBC_percent_MV1 <- sperm102s_mv1$DHA_in_RBC_._MV1
sperm102s$DHA_in_RBC_percent_MV4 <- sperm102s_mv4$DHA_in_RBC_._MV4
sperm102s$Change_in_percent_DHA_in_RBC <- sperm102s$DHA_in_RBC_percent_MV4 - sperm102s$DHA_in_RBC_percent_MV1

sperm102s %>% head()
sperm102s$EPA_in_Seminal_Fluid_MV1 <- sperm102s_mv1$EPA_in_Seminal_Fluid_._MV1
sperm102s$EPA_in_Seminal_Fluid_MV4 <- sperm102s_mv4$EPA_in_Seminal_Fluid_._MV4
sperm102s$Change_in_percent_EPA_in_Seminal_Fluid <- sperm102s$EPA_in_Seminal_Fluid_MV4 - sperm102s$EPA_in_Seminal_Fluid_MV1

sperm102s %>% head()
sperm102s$DHA_in_Seminal_Fluid_MV1 <- sperm102s_mv1$DHA_in_Seminal_Fluid_._MV1
sperm102s$DHA_in_Seminal_Fluid_MV4 <- sperm102s_mv4$DHA_in_Seminal_Fluid_._MV4
sperm102s$Change_in_percent_DHA_in_Seminal_Fluid <- sperm102s$DHA_in_Seminal_Fluid_MV4 - sperm102s$DHA_in_Seminal_Fluid_MV1

sperm102s %>% head()
sperm102s$EPA_in_Sperm_MV1 <- sperm102s_mv1$EPA_in_Sperm._MV1
sperm102s$EPA_in_Sperm_MV4 <- sperm102s_mv4$EPA_in_Sperm._MV4
sperm102s$Change_in_percent_EPA_in_Sperm <- sperm102s$EPA_in_Sperm_MV4 - sperm102s$EPA_in_Sperm_MV1

sperm102s %>% head()
sperm102s$DHA_in_Sperm_MV1 <- sperm102s_mv1$DHA_in_Sperm_._MV1
sperm102s$DHA_in_Sperm_MV4 <- sperm102s_mv4$DHA_in_Sperm_._MV4
sperm102s$Change_in_percent_DHA_in_Sperm <- sperm102s$DHA_in_Sperm_MV4 - sperm102s$DHA_in_Sperm_MV1

sperm102s$Group <- gsub("placebo", "Control", sperm102s$Group)
sperm102s$Group <- gsub("treatment", "Intervention", sperm102s$Group)

# BOXPLOT (CHANGE IN CONC OF VITD CONC SERUM)

stat.test.17s <- sperm17s %>%
  t_test(Change_in_VitD_Conc_Serum_nmol_L ~ Group) 
stat.test.17s 

stat.test.102s <- sperm102s %>%
  t_test(Change_in_VitD_Conc_Serum_nmol_L ~ Group) 
stat.test.102s

a <- ggplot(sperm17s, aes(x = Group, y = Change_in_VitD_Conc_Serum_nmol_L)) + 
      geom_boxplot() + 
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
      ylab("") +
      xlab("") +
      ylim(c(-50,250)) + 
      ggtitle("(ii) N = 17") +
      stat_pvalue_manual(stat.test.17s, y.position = 200,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 1, size = 5) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

b <- ggplot(sperm102s, aes(x = Group, y = Change_in_VitD_Conc_Serum_nmol_L)) + 
      geom_boxplot() + 
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
      ylab("Change in Conc. of Vitamin D \n in Serum (nmol/L)") +
      xlab("") + 
      ylim(c(-50,250)) + 
      ggtitle("(i) N = 102") +
      stat_pvalue_manual(stat.test.102s, y.position = 200,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 1, size = 5) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5) ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

ggarrange(b, a, ncol = 2, nrow = 1) %>% 
  ggexport(filename = paste0(dir,"/output/Figure3A_change_in_conc_vitD_in_serum_nmol_L.pdf"), width = 10,
           height = 5)

# BOXPLOT (CHANGE IN PERCENT OF EPA IN RBC)

stat.test.17s <- sperm17s %>%
  t_test(Change_in_percent_EPA_in_RBC ~ Group) 
stat.test.17s 

stat.test.102s <- sperm102s %>%
  t_test(Change_in_percent_EPA_in_RBC ~ Group) 
stat.test.102s

a <- ggplot(sperm17s, aes(x = Group, y = Change_in_percent_EPA_in_RBC)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("") +
  xlab("") +
  ylim(c(-2,3)) +
  ggtitle("(ii) N = 17") + 
  stat_pvalue_manual(stat.test.17s, y.position = 2,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.2, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

b <- ggplot(sperm102s, aes(x = Group, y = Change_in_percent_EPA_in_RBC)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("Change in % EPA in RBC \n % of total fatty acids") +
  xlab("") + 
  ylim(c(-2,3)) +
  ggtitle("(i) N = 102") + 
  stat_pvalue_manual(stat.test.102s, y.position = 2,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.5, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5) ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

ggarrange(b, a, ncol = 2, nrow = 1) %>% 
  ggexport(filename = paste0(dir,"/output/Figure3B_change_in_percent_EPA_in_RBC.pdf"), width = 10,
           height = 5)

# BOXPLOT (CHANGE IN PERCENT OF DHA IN RBC)

stat.test.17s <- sperm17s %>%
  t_test(Change_in_percent_DHA_in_RBC ~ Group) 
stat.test.17s 

stat.test.102s <- sperm102s %>%
  t_test(Change_in_percent_DHA_in_RBC ~ Group) 
stat.test.102s

a <- ggplot(sperm17s, aes(x = Group, y = Change_in_percent_DHA_in_RBC)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("") +
  xlab("") +
  ylim(c(-2,5.8)) +
  ggtitle("(ii) N = 17") + 
  stat_pvalue_manual(stat.test.17s, y.position = 4.5,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.5, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

b <- ggplot(sperm102s, aes(x = Group, y = Change_in_percent_DHA_in_RBC)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("Change in % DHA in RBC \n % of total fatty acids") +
  xlab("") + 
  ylim(c(-2,5.8)) +
  ggtitle("(i) N = 102") + 
  stat_pvalue_manual(stat.test.102s, y.position = 4.5,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.5, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5) ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

ggarrange(b, a, ncol = 2, nrow = 1) %>% 
  ggexport(filename = paste0(dir,"/output/Figure3C_change_in_percent_DHA_in_RBC.pdf"), width = 10,
           height = 5)
		   
# BOXPLOT (CHANGE IN PERCENT OF EPA IN SEMINAL FLUID)

stat.test.17s <- sperm17s %>%
  t_test(Change_in_percent_EPA_in_Seminal_Fluid ~ Group) 
stat.test.17s 

stat.test.102s <- sperm102s %>%
  t_test(Change_in_percent_EPA_in_Seminal_Fluid ~ Group) 
stat.test.102s

a <- ggplot(sperm17s, aes(x = Group, y = Change_in_percent_EPA_in_Seminal_Fluid)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("") +
  xlab("") +
  ylim(c(-0.5,0.5)) +
  ggtitle("(ii) N = 17") + 
  stat_pvalue_manual(stat.test.17s, y.position = 0.3,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.05, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

b <- ggplot(sperm102s, aes(x = Group, y = Change_in_percent_EPA_in_Seminal_Fluid)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("Change in % EPA in Seminal Fluid \n % of total fatty acids") +
  xlab("") + 
  ylim(c(-0.5,0.5)) +
  ggtitle("(i) N = 102") + 
  stat_pvalue_manual(stat.test.102s, y.position = 0.3,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.05, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5) ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

ggarrange(b, a, ncol = 2, nrow = 1) %>% 
  ggexport(filename = paste0(dir,"/output/Figure3D_change_in_percent_EPA_in_Seminal_fluid.pdf"), width = 10,
           height = 5)

# BOXPLOT (CHANGE IN PERCENT OF EPA IN SPERM)

stat.test.17s <- sperm17s %>%
  t_test(Change_in_percent_EPA_in_Sperm ~ Group) 
stat.test.17s 

stat.test.102s <- sperm102s %>%
  t_test(Change_in_percent_EPA_in_Sperm ~ Group) 
stat.test.102s

a <- ggplot(sperm17s, aes(x = Group, y = Change_in_percent_EPA_in_Sperm)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("") +
  xlab("") +
  ylim(c(-0.25,0.25)) +
  ggtitle("(ii) N = 17") + 
  stat_pvalue_manual(stat.test.17s, y.position = 0.1,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.01, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

b <- ggplot(sperm102s, aes(x = Group, y = Change_in_percent_EPA_in_Sperm)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  ylab("Change in % EPA in Sperm \n % of total fatty acids") +
  xlab("") + 
  ylim(c(-0.25,0.25)) +
  ggtitle("(i) N = 102") + 
  stat_pvalue_manual(stat.test.102s, y.position = 0.1,label = "{format(p, digits = 2, scientific = T)}",vjust = -1, bracket.nudge.y = 0.01, size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5) ,text = element_text(size = 15,color="black"), axis.text =element_text(color="black", size = 15))

ggarrange(b, a, ncol = 2, nrow = 1) %>% 
  ggexport(filename = paste0(dir,"/output/Figure3E_change_in_percent_EPA_in_Sperm.pdf"), width = 10,
           height = 5)
		   		   
