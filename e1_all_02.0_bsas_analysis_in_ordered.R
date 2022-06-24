## set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## LOAD LIBRARIES ----

library(betareg) 
library(car) 
require(devtools) 
library(grid)
library(gridExtra)
library(ggpubr)
library(vegan)
library(data.table)
library(broom)
library(viridis)
library(multcomp)
library(tidyverse)

## IMPORT AND FORMAT DATA ----

## import data ----

## import raw data from geneious variant calling (3 miseq runs, 2 time points)
## for this experiment 
data_6mo_run1 <- read.csv("in/c3_6.1_02.1_geneious_data.csv") 
data_6mo_run2 <- read.csv("in/c3_6.2_02.1_geneious_data.csv")
data_12mo <- read.csv("in/c3_12_02.1_geneious_data.csv")

## import fish data
c3_fish_data <- read.csv("in/c3_all_02.3_fish_data.csv") 

## import amplicon file
c3_amplicon_cpg <- read.csv("in/c3_all_02.2_gene_cpg_dmrt1_split.csv")

## import Geneious data from Budd et al., 2022 (for adult data)
c2_data_unmod <- read.csv("in/c2_02.1_geneious_data_2021.csv")

## import sampling info (location etc.) from Budd et al., 2022
c2_field_data <- read.csv("in/c2_02.4_field_data.csv")

## add fish ID numbers from Budd et al., 2022
c2_fish_ids=read.csv("in/c2_02.3_fish_ids.csv")

## import Geneious data from Domingos et al., 2019 (for adult data)
domingos_data <- read.csv("in/domingos_et_al_2018_data.csv")

## format data from this experiment (aka chapter 3) ----

## combine 6 and 12 month data for experiment 1
c3_data <- rbind(data_6mo_run1, data_6mo_run2, data_12mo) %>%
  ## remove annotations that are not polymorphisms
  filter(Type=="Polymorphism") %>% 
  ## remove non 'C-T' changes 
  filter(Change=="C -> T") %>%
  ## add amplicon and CpG names
  mutate(Minimum = as.numeric(Minimum)) %>%
  left_join(., c3_amplicon_cpg, by = "Minimum") %>%
  ## rename 'Position' column 
  rename(CpG_site = "Position") %>%
  ## add fish numbers, treatments, lengths and weights
  left_join(., c3_fish_data %>%
              select(FishNo, Treatment, Month_Treatment, 
                     Month, Length, Weight, TrackName) %>%
              rename(Track.Name = TrackName), 
            by = "Track.Name") %>%
  ## remove NAs
  na.omit(CpG_site) %>% 
  droplevels(.) %>%
  ## remove % signs and create proportion column
  mutate(Reference.Frequency = sub("\\%.*","\\", Reference.Frequency)) %>%
  mutate(Reference.Frequency = as.numeric(as.character(Reference.Frequency))) %>%
  mutate(prop_meth = Reference.Frequency/100) %>%
  ## delete unwanted columns
  select(-c(Variant.Raw.Frequency, Variant.Frequency, 
            Name, Change, Type, Track.Name, 
            Minimum, Maximum, Coverage, Reference.Frequency,
            Length, Weight))

## format adult data from Budd et al., 2022 (aka chapter 2) ----

c2_data <- c2_data_unmod %>% 
  ## drop gene annotation information from 'type' column
  filter(Type!="gene") %>% 
  ## remove annotations that are not polymorphisms
  filter(Type=="Polymorphism") %>% 
  ## remove SNPs that didn't occur at CpG sites
  na.omit(CpG_site) %>%
  ## remove non-C-T changes
  filter(Change=="C -> T") %>%
  ## remove anything less than 500 reads 
  filter(., Coverage>499) %>%
  ## add gene and cpg information
  mutate(Minimum = as.numeric(Minimum)) %>%
  left_join(., c3_amplicon_cpg, by = "Minimum") %>%
  ## add fish IDs 
  left_join(., c2_fish_ids, by = "Track.Name") %>%
  rename(Fish = LTMPno) %>%
  ## remove % signs and calculate prop meth
  mutate(Reference.Frequency = sub("\\%.*","\\", Reference.Frequency)) %>%
  mutate(Reference.Frequency = as.numeric(as.character(Reference.Frequency))) %>%
  mutate(prop_meth = Reference.Frequency/100) %>%
  ## remove barra from low frequency catch locations 
  filter(., Fish!="#251") %>%
  filter(., Fish!="#254") %>%
  filter(., Fish!="#110") %>%
  ## add sampling location info
  left_join(., c2_field_data %>% 
              select(JulieRegion, SexCode, Sex, tl, 
                     yr, AgeClass, YearClass, MonthCaught, LTMPno) %>%
              rename(Fish = LTMPno), 
            by = "Fish") %>%
  ## remove intersex fish (n = 1)
  filter(., Sex!="Transitional") %>% 
  droplevels(.) %>%
  ## delete unwanted columns
  select(-c(Variant.Raw.Frequency, Variant.Frequency, 
            Name, Change, Type, Track.Name)) %>%
  ## rename some others
  rename(Total_length = tl, 
         Percent_methylation = Reference.Frequency) %>%
  ## create size bins
  mutate(Length_binned = cut(Total_length, 
                             breaks = c(0,60,70,80,90,100,120), 
                             labels = c("50-60","60-70","70-80","80-90","90-100","v100"))) %>%
  ## make factors factors
  mutate(Sex = as.factor(Sex), 
         Gene = as.factor(Gene), 
         JulieRegion = as.factor(JulieRegion), 
         CpG_site = Position) %>%
  ## edit region and sex names so that east coast males come first alphabetically
  mutate(JulieRegionEC = gsub("Mid-northern GoC", "MidNorthernGulf", JulieRegion) %>%
           gsub("Southern GoC", "SouthernGulf", .) %>%
           gsub("Wet-tropics East Coast", "EastCoastWetTropics", .)) %>%
  mutate(SexM = gsub("Female", "xFemale", Sex)) %>%
  ## alternatively, just reorder the factor here (although df edits can cause headaches)
  mutate(JulieRegionEC = as.factor(JulieRegionEC)) %>%
  mutate(SexM = as.factor(SexM))

## use east coast barra only
c2_data_ec<-subset(c2_data, JulieRegionEC=="EastCoastWetTropics")

## summarise age for east coast males
age_ltmp_ec_males <- unique(select(c2_data_ec, Fish, Sex, AgeClass)) %>%
  filter(., Sex == "Male") %>%
  summarise(., mean = mean(AgeClass), SD = sd(AgeClass),
            min = min(AgeClass), max = max(AgeClass))

## summarise age for east coast females 
age_ltmp_ec_females <- unique(select(c2_data_ec, Fish, Sex, AgeClass)) %>%
  filter(., Sex == "Female")  %>%
  summarise(., mean = mean(AgeClass), SD = sd(AgeClass),
            min = min(AgeClass), max = max(AgeClass))

deletes<-c("Minimum", "Maximum", "Coverage", "Percent_methylation",
           "JulieRegion", "yr", "AgeClass", "YearClass", 
           "MonthCaught", "JulieRegionEC", "SexCode", "SexM")

c2_data_ec <- c2_data_ec[ , !(names(c2_data_ec) %in% deletes)]

## add Gene column
c2_data_ec$Gene <- c2_data_ec$Amplicon
c2_data_ec$Gene <- gsub("esr1_PE1", "esr1", c2_data_ec$Gene)
c2_data_ec$Gene <- gsub("cyp19a1a_PE1", "cyp19a1a", c2_data_ec$Gene)
c2_data_ec$Gene <- gsub("nr5a2_PE1", "nr5a2", c2_data_ec$Gene)
c2_data_ec$Gene <- gsub("dmrt1_PE1_C", "dmrt1", c2_data_ec$Gene) 

## edit names
names(c2_data_ec)[names(c2_data_ec)=="Fish"]<-"FishNo" 
names(c2_data_ec)[names(c2_data_ec)=="Sex"]<-"Treatment" 

## edit month and treatment info
c2_data_ec$Month<-"Wild"
c2_data_ec$Month<-as.factor(c2_data_ec$Month)
c2_data_ec$Month_Treatment <- paste(c2_data_ec$Treatment, c2_data_ec$Month, sep=".")

## retain only relevant columns
c2_data_ec <- c2_data_ec[c("Amplicon", "CpG_site", "FishNo", "Treatment", 
                           "Month_Treatment", "Month", "prop_meth", "Gene")]

## make month a factor for c3, too
c3_data$Month<-as.factor(c3_data$Month)

## combine data sets
data_both <- rbind.data.frame(c2_data_ec, c3_data)

## it seems that there is a lot of unecessary factor-making going on in this script...
data_both$Amplicon <- factor(data_both$Amplicon)
data_both$Month<-factor(data_both$Month, levels=c("Wild", "12", "6"))

## remove Gene column
data_both <- select(data_both, -Gene)

## format adult data from Domingos et al., 2018 ----

## keep wild fish only
domingos_data_wild<-filter(domingos_data, Month == "Wild") %>%
  rename(., prop_meth = Reference.Frequency.P)

data_all<-rbind.data.frame(domingos_data_wild, data_both) 
data_all$Amplicon <- factor(data_all$Amplicon)

data_all$Month<-factor(data_all$Month, 
                       levels=c("Wild", "12", "6"))  

## format combined data frame ----

## rename adult fish and amplicons 
data_all$Month <- gsub("Wild", "Adult", data_all$Month)
data_all$Amplicon <- gsub("dmrt1_PE1_C", "dmrt1_PE1", data_all$Amplicon)
data_all$Amplicon <- gsub("dmrt1_Ppt1", "dmrt1_P", data_all$Amplicon)
data_all$Amplicon <- gsub("amh_E7o7", "amh_E7", data_all$Amplicon)
data_all$Amplicon <- gsub("cyp19a1a_E8o9", "cyp19a1a_E8", data_all$Amplicon)
data_all$Amplicon <- gsub("foxl2_c3pt1", "foxl2_c3.1", data_all$Amplicon)
data_all$Amplicon <- gsub("foxl2_c3pt2", "foxl2_c3.2", data_all$Amplicon)
data_all$Amplicon <- as.factor(data_all$Amplicon)

## change to HT LT etc 
data_all$Treatment <- gsub("T2434", "FT", data_all$Treatment)
data_all$Treatment <- gsub("T24", "LT", data_all$Treatment)
data_all$Treatment <- gsub("T34", "HT", data_all$Treatment)
data_all$Treatment <- gsub("T29", "Control", data_all$Treatment)

## edit month labels
data_all$Month <- gsub("6", "6 mph", data_all$Month)
data_all$Month <- gsub("12", "12 mph", data_all$Month)
data_all$Month <- gsub("Wild", "Adult", data_all$Month)

## CREATE FIGURES ----

## make new data frame for means (because ggplot2 plots stat identity)
data_all_sum <- data_all %>% 
  group_by(Amplicon, Month_Treatment, Month, Treatment, CpG_site) %>% 
  summarise(., mean_prop_meth = mean(prop_meth, na.rm=TRUE)) %>%
  ## convert back to percentage
  mutate(mean_perc_meth = mean_prop_meth*100)

## create df with means accross all (-1000; 'All' in figure)
data_cpg_means <- data_all_sum %>% 
  rename(avg_perc_meth = mean_perc_meth, avg_prop_meth = mean_prop_meth) %>%
  group_by(Amplicon, Month_Treatment, Month, Treatment) %>% 
  summarise(., mean_perc_meth = mean(avg_perc_meth, na.rm=TRUE),
            mean_prop_meth = mean(avg_prop_meth, na.rm=TRUE)) %>%
  mutate(CpG_site = -1000)

## merge the two
data_all_sum <- rbind.data.frame(data_all_sum, data_cpg_means)

## factors  
data_all_sum$Month <- factor(data_all_sum$Month, levels=c("6 mph", "12 mph", "Adult"))
data_all_sum$Month_Treatment <- as.factor(data_all_sum$Month_Treatment)

## logit transform data 
logitTransform <- function(p) { log(p/(1-p)) }
data_all$logitMeth <- logitTransform(data_all$prop_meth)

## main figure (PE1 heatmap grid) ----

PE1s<-droplevels(unique(data_all_sum$Amplicon[grep("*._PE1", data_all_sum$Amplicon)]))

for (a in levels(PE1s)) {
  data_all_amplicon <- filter(data_all_sum, Amplicon==a)
  data_all_amplicon$Treatment <- factor(data_all_amplicon$Treatment, 
                                        levels=c("Control", "LT", "FT", "HT", "Female", "Male"))
  data_temp <- filter(data_all, Amplicon==a)
  data_temp$Month_Treatment <- as.factor(data_temp$Month_Treatment)
  data_all_amplicon$CpG_site <- as.factor(data_all_amplicon$CpG_site) 
  name_temp = paste("HM", a, sep="_")
  print(name_temp)
  ## anova
  anova_temp <- aov(logitMeth ~ Month_Treatment, data = data_temp)
  glht_temp <- summary(glht(anova_temp, linfct = mcp(Month_Treatment = "Tukey")), 
                       test = adjusted("bonferroni"))
  cld<-cld(glht_temp)
  cld2<-as.data.frame(cld$mcletters$Letters, rownames=NULL)
  cld2<-setDT(cld2, keep.rownames = TRUE)[]
  names(cld2)[names(cld2)=="rn"]<-"Month_Treatment" 
  names(cld2)[names(cld2)=="cld$mcletters$Letters"]<-"Letters"
  data_all_amplicon <- merge(data_all_amplicon, cld2, by  = "Month_Treatment")
  ## add tukey letters
  data_all_amplicon$ylabs <- paste0(data_all_amplicon$Treatment, " (", 
                                   data_all_amplicon$Letters, ")")
  ## order by treatment
  ordered<-c(data_all_amplicon$ylabs[grep("Control.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("LT.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("FT.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("HT.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("Female.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("Male.*", data_all_amplicon$ylabs)])
  data_all_amplicon$ylabs <- factor(data_all_amplicon$ylabs, levels=unique(ordered))
  ## make x labs vector (to relabel -1000 as 'All')
  xlabs <- c("All", 
             levels(factor(data_all_amplicon$CpG_site))[levels(factor(data_all_amplicon$CpG_site))!="-1000"])
  ## make plot
  HM<-ggplot(data = data_all_amplicon, 
              aes(x=as.character(CpG_site), 
                  y=ylabs, 
                  fill=mean_perc_meth)) + 
    geom_tile()+
    ggtitle(gsub("_.*", "", a))+
    xlab("CpG position (bp)")+
    ylab("Treatment")+
    scale_x_discrete(labels=xlabs)+
    theme(axis.text.x = element_text(size=9, angle = 45, hjust = 1),
          axis.text.y = element_text(size=9),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          strip.text.y = element_text(size = 9, margin = margin(0, 0, 0, 0, "cm")),
          panel.border = element_blank(),
          strip.background = element_blank(), 
          plot.margin = unit(c(0.1,0.2,0.1,0.2), "cm"), ## t, r, b, l
          legend.margin = margin(0,0,0,0),
          legend.title = element_text(size = 10),
          legend.key.width = unit(0.5, "line"),
          legend.text = element_text(size=8),
          plot.title = element_text(size=11, face = "bold.italic", hjust=0))+ 
    geom_vline(xintercept=1.5, colour="grey10", size=0.3)+
    # guides(fill = guide_legend(override.aes = list(size = 0.5)))
  facet_grid(rows=vars(Month), scales = "free", space = "free_y")
  HM_l <- HM + scale_fill_viridis(direction = -1, discrete=FALSE, option="G", 
                                      limits=c(-10, 100), 
                                      name = "DNAm\n(%)")
  HM_nl <- HM + scale_fill_viridis(direction = -1, discrete=FALSE, option="G",
                                  name = "DNAm\n(%)")
  assign(name_temp, HM_l)
  assign(paste0(name_temp, "_nl"), HM_nl)
}

## combine no limits plot (no shared legend)
figure_main_nl <- ggarrange(HM_nr5a2_PE1_nl+theme(axis.title.x = element_blank()),
                      HM_cyp19a1a_PE1_nl+theme(axis.title.y = element_blank(), 
                                               axis.title.x = element_blank()),
                      HM_dmrt1_PE1_nl+theme(axis.title.x = element_blank()),
                      HM_esr1_PE1_nl+theme(axis.title.y = element_blank(), 
                                           axis.title.x = element_blank()),
                      HM_sox9_PE1_nl,
                      HM_amh_PE1_nl+theme(axis.title.y = element_blank()),
                      heights = c(rep(c(0.96, 0.96, 1), 2)),
                      nrow=3, ncol=2)

figure_main_nl

## save figure (no limits - i.e., not scaled)
ggsave("out/figure_main_nl.pdf", figure_main_nl, 
       width=21*1.1, height=29.7*0.75, units="cm", device="pdf")

## combine WITH limits plot (no shared legend)
figure_main_l <- ggarrange(HM_nr5a2_PE1+theme(axis.title.x = element_blank()),
                         HM_cyp19a1a_PE1+theme(axis.title.y = element_blank(), 
                                               axis.title.x = element_blank()),
                         HM_dmrt1_PE1+theme(axis.title.x = element_blank()),
                         HM_esr1_PE1+theme(axis.title.y = element_blank(), 
                                           axis.title.x = element_blank()),
                         HM_sox9_PE1,
                         HM_amh_PE1+theme(axis.title.y = element_blank()),
                         heights = c(rep(c(0.96, 0.96, 1), 2)),
                         nrow=3, ncol=2
)

figure_main_l

## save figure (WITH limits - i.e., scaled 0 - 100)
ggsave("out/figure_main_l.pdf", figure_main_l, 
       width=21*1.1, height=29.7*0.75, units="cm", device="pdf")

## supp figures (individual heatmaps) ----

for (a in levels(factor(data_all$Amplicon))) {
  data_all_amplicon <- filter(data_all_sum, Amplicon==a)
  data_temp <- filter(data_all, Amplicon==a)
  data_temp$Month_Treatment <- as.factor(data_temp$Month_Treatment)
  data_all_amplicon$CpG_site <- factor(data_all_amplicon$CpG_site) 
  name_temp = paste("HM", a, sep="_")
  print(name_temp)
  ## anova
  anova_temp <- aov(logitMeth ~ Month_Treatment, data = data_temp)
  glht_temp <- summary(glht(anova_temp, linfct = mcp(Month_Treatment = "Tukey")), 
                       test = adjusted("bonferroni"))
  cld<-cld(glht_temp)
  cld2<-as.data.frame(cld$mcletters$Letters, rownames=NULL)
  cld2<-setDT(cld2, keep.rownames = TRUE)[]
  names(cld2)[names(cld2)=="rn"]<-"Month_Treatment" 
  names(cld2)[names(cld2)=="cld$mcletters$Letters"]<-"Letters"
  data_all_amplicon <- merge(data_all_amplicon, cld2, by  = "Month_Treatment")
  ## add tukey letters
  data_all_amplicon$ylabs <- paste0(data_all_amplicon$Treatment, " (", 
                                    data_all_amplicon$Letters, ")")
  ## order by treatment
  ordered<-c(data_all_amplicon$ylabs[grep("Control.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("LT.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("FT.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("HT.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("Female.*", data_all_amplicon$ylabs)], 
             data_all_amplicon$ylabs[grep("Male.*", data_all_amplicon$ylabs)])
  data_all_amplicon$ylabs <- factor(data_all_amplicon$ylabs, levels=unique(ordered))
  ## make x labs vector (to relabel -1000 as 'All')
  xlabs <- c("All", 
             levels(factor(data_all_amplicon$CpG_site))[levels(factor(data_all_amplicon$CpG_site))!="-1000"])
  ## make plot
  HM<-ggplot(data = data_all_amplicon, 
             aes(x=as.character(CpG_site), 
                 y=ylabs, 
                 fill=mean_perc_meth)) + 
    geom_tile()+
    ggtitle(gsub("_", " ", a))+
    xlab("CpG position (bp)")+
    ylab("Treatment")+
    scale_x_discrete(labels = xlabs)+
    theme(axis.text.x = element_text(size=9, angle = 45, hjust = 1),
          axis.text.y = element_text(size=9),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          strip.text.y = element_text(size = 9, margin = margin(0, 0, 0, 0, "cm")),
          panel.border = element_blank(),
          strip.background = element_blank(), 
          plot.margin = unit(c(0.1,0.2,0.1,0.2), "cm"),
          legend.margin = margin(0,0,0,0),
          legend.title = element_text(size = 10),
          legend.key.width = unit(0.5, "line"),
          legend.text = element_text(size=8),
          plot.title = element_text(size=11, face = "bold.italic", hjust=0))+ 
    theme(legend.position = "none")+ #remove legend
    geom_vline(xintercept=1.5, colour="grey10", size=0.3)+
    facet_grid(rows=vars(Month), scales = "free", space = "free_y")
  HM_l <- HM + scale_fill_viridis(direction = -1, discrete=FALSE, option="G", 
                                  limits=c(0, 110), 
                                  name = "DNAm\n(%)")
  ggsave(plot=HM_l, filename=paste0("out/", name_temp, "_l.pdf"), 
         device="pdf", width=16, height=11, units=c("cm"))
}

## OTHER STATISTICS ----

## create mean +/- SD, min, max table ----
data_all_msd = data_all %>% 
  group_by(Amplicon, Month_Treatment, Month, Treatment) %>% 
  summarise(mean_RFP = mean(prop_meth, na.rm=TRUE),
            sd_RFP = sd(prop_meth, na.rm=TRUE), 
            min_RFP = min(prop_meth, na.rm=TRUE), 
            max_RFP = max(prop_meth, na.rm=TRUE)) %>%
  mutate('Mean (+/- SD)' = paste0(round(mean_RFP*100, 2), 
                                  " (", SD = round(sd_RFP*100, 2), ")"),
         Min. = round(min_RFP*100, 2), 
         Max. = round(max_RFP*100, 2), 
         Amplicon = gsub("_", " ", Amplicon)) %>%
  ungroup(.) %>%
  select(-c(mean_RFP, min_RFP, max_RFP, sd_RFP, Month_Treatment)) %>%
  rename('Fish age' = Month)

## export it
write.csv(data_all_msd, "out/c3_meth_summary_table.csv", row.names=FALSE)

## create glht table ----

glht_table<-NULL

for (a in levels(data_all_sum$Amplicon)) {
  data_all_amplicon <- filter(data_all_sum, Amplicon==a)
  data_temp <- filter(data_all, Amplicon==a)
  data_temp$Month_Treatment <- as.factor(data_temp$Month_Treatment)
  data_all_amplicon$CpG_site <- factor(data_all_amplicon$CpG_site) 
  print(a)
  ## anova
  anova_temp <- aov(logitMeth ~ Month_Treatment, data = data_temp)
  glht_temp <- summary(glht(anova_temp, linfct = mcp(Month_Treatment = "Tukey")), 
                       test = adjusted("bonferroni"))
  cld<-cld(glht_temp)
  ## table
  glht_temp2<-tidy(summary(glht_temp))
  glht_temp3<-cbind(a, glht_temp2)
  glht_table<-rbind(glht_table, glht_temp3)
}

## edit glht table
glht_table_edit <- glht_table %>%
  mutate(contrast = gsub("T2434", "FT", contrast) %>%
           gsub("T24", "LT", .) %>%
           gsub("T34", "HT", .) %>%
           gsub("T29", "Control", .) %>%
           gsub("_", " mph ", .) %>%
           gsub("Female.Wild", "Adult Female", .) %>%
           gsub("Male.Wild", "Adult Male", .)) %>%
  mutate(a = gsub("_", " ", a), 
         estimate = round(estimate, digits = 3), 
         std.error = round(std.error, digits = 3), 
         statistic = round(statistic, digits = 3), 
         adj.p.value = round(adj.p.value, digits = 3)) %>%
  select(-c(term, null.value)) %>%
  mutate('Sig. Code' = ifelse(adj.p.value <=0.001, "***",
                              ifelse(adj.p.value<=0.01, "**",
                                     ifelse(adj.p.value <=0.05, "*",
                                            ifelse(adj.p.value <=0.1, ".",
                                                   " "))))) %>%
  rename(Comparison = contrast, Amplicon = a, 
         Estimate = estimate, 'Std. error' = std.error, 
         Statistic = statistic, 'Adj. p-value' = adj.p.value)

## export glht table
write.csv(glht_table_edit, "out/c3_glht_table.csv", row.names=FALSE)

## within gene comparisons (inc. exon t-tests) ----

data_all_less <- data_all %>% 
  ## remove adults (no 3' exon data)
  filter(., Month != "Adult") %>%
  ## convert to percentage
  mutate(perc_meth = prop_meth*100) 

for (i in levels(factor(data_all_less$Amplicon))) {
  # print(paste0("data_", i))
  assign(paste0("data_", i), filter(data_all_less, Amplicon == i))
}

max(data_foxl2_P$perc_meth, na.rm=TRUE)
max(data_foxl2_E1pt1$perc_meth, na.rm=TRUE)
max(data_foxl2_E1pt2$perc_meth, na.rm=TRUE)

mean(data_dmrt1_P$perc_meth, na.rm=TRUE)
mean(data_dmrt1_PE1$perc_meth, na.rm=TRUE)
mean(data_dmrt1_E5$perc_meth, na.rm=TRUE)

mean(data_amh_E7$perc_meth, na.rm=TRUE)
mean(data_cyp19a1a_E8$perc_meth, na.rm=TRUE)
mean(data_dmrt1_E5$perc_meth, na.rm=TRUE)

amh<-t.test(data_amh_PE1$perc_meth, data_amh_E7$perc_meth, var.equal=FALSE) 
cyp19a1a<-t.test(data_cyp19a1a_PE1$perc_meth, data_cyp19a1a_E8$perc_meth, var.equal=FALSE) 
dmrt1<-t.test(data_dmrt1_PE1$perc_meth, data_dmrt1_E5$perc_meth, var.equal=FALSE) 

paste0("amh: t(", round(amh$parameter[[1]], digits=2), 
       ") = ", round(amh$statistic[[1]], digits=2), 
       ", p = ", ifelse(amh$p.value < 0.001, 
                      "<0.001", round(amh$p.value, digits=2)))

paste0("cyp19a1a: t(", round(cyp19a1a$parameter[[1]], digits=2), 
       ") = ", round(cyp19a1a$statistic[[1]], digits=2), 
       ", p = ", ifelse(cyp19a1a$p.value < 0.001, 
                      "<0.001", round(cyp19a1a$p.value, digits=2)))

paste0("dmrt1: t(", round(dmrt1$parameter[[1]], digits=2), 
       ") = ", round(dmrt1$statistic[[1]], digits=2), 
       ", p = ", ifelse(dmrt1$p.value < 0.001, 
                      "<0.001", round(dmrt1$p.value, digits=2)))

## end script 

