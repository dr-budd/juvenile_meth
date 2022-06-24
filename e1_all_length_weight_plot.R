## set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(multcomp)
library(data.table) ## setDT
library(viridis)
library(car)     
library(tidyverse)

## IMPORT AND FORMAT DATA ----

## import data
length_weight_data_unmod <- read.csv("in/c3_length_weight_data_all.csv")

## format data
length_weight_data <- length_weight_data_unmod %>%
  mutate(Treatment = gsub("F", "FT", Treatment) %>%
           gsub("24", "LT", .) %>%
           gsub("34", "HT", .) %>%
           gsub("29", "Control", .)) %>%
  mutate(Treatment = factor(Treatment, levels=c("Control", "LT", "FT", "HT")), 
         Months = as.factor(Months)) %>%
  mutate(Month = paste0(Months, " mph")) %>%
  mutate(Month = factor(Month, levels=c("12 mph", "9 mph", "6 mph")))

## WEIGHT STATS ----

## se function
se <- function(x) sqrt(var(x)/length(x))

## calculate mean SD/SE
length_weight_data_msdse = length_weight_data %>% 
  group_by(Months, Treatment) %>% 
  summarise(mean_len = mean(Length, na.rm=T),
            sd_len = sd(Length, na.rm=T), 
            mean_wei = mean(Weight, na.rm=T), 
            sd_wei = sd(Weight, na.rm=T), 
            mean_len = mean(Length, na.rm=T),
            se_len = se(Length), 
            mean_wei = mean(Weight, na.rm=T), 
            se_wei = se(Weight)) %>%
            ungroup(.)%>%
            # mutate(prev_mean_wei = c(1, 1, 1, 1, length_weight_data_msdse$mean_wei[1:8])) %>%
            mutate(prev_mean_wei = c(0.01, 0.01, 0.01, 0.01, pull(., var=mean_wei)[1:8])) %>%
            mutate(days = c(rep(180, 4), rep(90, 8))) %>%
            mutate(SGR = (log(mean_wei) - log(prev_mean_wei))/days*100)

## make a nicely formatted supplementary table
## (did not end up using SGR)
w_supp_table <- length_weight_data_msdse %>%
  mutate("Mean length (cm)" = round(mean_len, 2)) %>%
  mutate("SD length" = round(sd_len, 2)) %>%
  mutate("SE length" = round(se_len, 2)) %>%
  mutate("Mean weight (g)" = round(mean_wei, 2)) %>%
  mutate("SD weight" = round(sd_wei, 2)) %>%
  mutate("SE weight" = round(se_wei, 2)) %>%
  # mutate("SGR (%)" = round(SGR, 2)) %>%
  select(-c("mean_len", "sd_len", "se_len", 
            "mean_wei", "sd_wei", "se_wei", 
            "prev_mean_wei", "SGR", "days"))

## export it
write.csv(w_supp_table, "out/c3_w_supp_table_both.csv", row.names = FALSE)

## run among treatment comparisons ----

## create empty df
cld_table_weight_369 <- NULL

## create tukey letters df
for (i in levels(length_weight_data$Month)) {
  temp_data<-filter(length_weight_data, Month==i)
  temp_name<-paste("Weight", i, "Months", sep="_")
  print(temp_name)
  temp_aov_w <- aov(Weight ~ Treatment, data = temp_data)
  temp_glht_w<-(glht(temp_aov_w, linfct = mcp(Treatment = "Tukey")))
  cld_w<-cld(temp_glht_w)
  cld_table_w_369<-as.data.frame(cld_w$mcletters$Letters, rownames=NULL)
  cld_table_w_369<-setDT(cld_table_w_369, keep.rownames = TRUE)[]
  cld_table_w_369$Month <- paste(i)
  names(cld_table_w_369)[names(cld_table_w_369)=="rn"]<-"Treatment" 
  names(cld_table_w_369)[names(cld_table_w_369)=="cld$mcletters$Letters"]<-"Letters" 
  letters_table_w_369 <- cld_table_w_369[match(paste(temp_data$Treatment),
                                       paste(cld_table_w_369$Treatment)),]
  ylabelpos <- NULL
  for (j in temp_data$Treatment) { 
    subsetis <- (temp_data[which(temp_data$Treatment==j),])
    ylabelpos <- c(ylabelpos,max(subsetis$Weight))
  }
  letters_table_w_369$ylabelpos <- ylabelpos*1.2
  cld_table_weight_369<-rbind(cld_table_weight_369, letters_table_w_369)
}

## edit table
cld_table_weight_369_unique<-(unique(cld_table_weight_369)) %>%
  rename(weight_label_letter = `cld_w$mcletters$Letters`, 
         weight_ylabelpos = ylabelpos) %>%
  mutate(Months = gsub(" mph", "", Month))

## add cld to main df
length_weight_data_msdse_tukey <- length_weight_data_msdse %>%
  left_join(., cld_table_weight_369_unique, by = c("Months", "Treatment")) %>%
  mutate(Months = factor(Months, levels=c("6", "9", "12"))) %>%
  mutate(Treatment = factor(Treatment, levels=c("Control", "LT", "FT", "HT"))) %>%
  ## order = control, LT, FT, HT
  mutate(weight_ylabelpos2 = c(c(200, 133, 100, 166), c(433, 466, 400, 500), c(566, 600, 633, 666))) %>%
  mutate(weight_xlabelpos = c(rep(0.8, 4), rep(1.6, 4), rep(3.2, 4)))

## CREATE PLOT ----

## edit viridis because yellow is hard to see
viridis_colours <- viridis(n = 4)
viridis_colours_darkyellow <- gsub("#FDE725FF", "darkgoldenrod2", viridis_colours)
viridis_colours_darkyellow

## plot weight with tukey letters
line_weight<-ggplot(length_weight_data_msdse_tukey, 
                    aes(factor(Months), mean_wei, colour=Treatment))+
  geom_line(aes(group=Treatment))+
  geom_point()+
  scale_colour_manual(values=viridis_colours_darkyellow)+
  geom_text(aes(weight_xlabelpos, weight_ylabelpos2, 
                 label = weight_label_letter), 
             fontface = "bold", size = 4, 
            show.legend = FALSE)+
  geom_errorbar(aes(ymin=mean_wei-se_wei, ymax=mean_wei+se_wei, group=Treatment),
                width=0.05)+
  theme_bw()+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.background = element_blank())+
  ylab("Weight (g)") +
  xlab("Time (mph)")

line_weight

## LENGTH STATS ----

## create empty df
cld_table_length_369 <- NULL

## create tukey letters df
for (i in levels(length_weight_data$Month)) {
  temp_data<-filter(length_weight_data, Month==i)
  temp_name<-paste("Length", i, "Months", sep="_")
  print(temp_name)
  temp_aov_w <- aov(Length ~ Treatment, data = temp_data)
  temp_glht_w<-(glht(temp_aov_w, linfct = mcp(Treatment = "Tukey")))
  cld_w<-cld(temp_glht_w)
  cld_table_w_369<-as.data.frame(cld_w$mcletters$Letters, rownames=NULL)
  cld_table_w_369<-setDT(cld_table_w_369, keep.rownames = TRUE)[]
  cld_table_w_369$Month <- paste(i)
  names(cld_table_w_369)[names(cld_table_w_369)=="rn"]<-"Treatment" 
  names(cld_table_w_369)[names(cld_table_w_369)=="cld$mcletters$Letters"]<-"Letters" 
  letters_table_w_369 <- cld_table_w_369[match(paste(temp_data$Treatment),
                                               paste(cld_table_w_369$Treatment)),]
  ylabelpos <- NULL
  for (j in temp_data$Treatment) { 
    subsetis <- (temp_data[which(temp_data$Treatment==j),])
    ylabelpos <- c(ylabelpos,max(subsetis$Length))
  }
  letters_table_w_369$ylabelpos <- ylabelpos*1.2
  cld_table_length_369<-rbind(cld_table_length_369, letters_table_w_369)
}

## edit table
cld_table_length_369_unique<-(unique(cld_table_length_369)) %>%
  rename(length_label_letter = `cld_w$mcletters$Letters`, 
         length_ylabelpos = ylabelpos) %>%
  mutate(Months = gsub(" mph", "", Month))

## add cld to main df
length_weight_data_msdse_tukey <- length_weight_data_msdse %>%
  left_join(., cld_table_length_369_unique, by = c("Months", "Treatment")) %>%
  mutate(Months = factor(Months, levels=c("6", "9", "12"))) %>%
  mutate(Treatment = factor(Treatment, levels=c("Control", "LT", "FT", "HT"))) %>%
  # order = control, LT, FT, HT
  mutate(length_ylabelpos2 = c(c(23, 21, 20, 22), c(29, 31, 30, 32), c(34, 35, 36, 37))) %>%
  mutate(length_xlabelpos = c(rep(0.8, 4), rep(1.6, 4), rep(3.2, 4))) 

## plot weight with tukey letters
line_length<-ggplot(length_weight_data_msdse_tukey, 
                    aes(factor(Months), mean_len, colour=Treatment))+
  geom_line(aes(group=Treatment))+
  geom_point()+
  scale_colour_manual(values=viridis_colours_darkyellow)+
  geom_text(aes(length_xlabelpos, length_ylabelpos2, 
                label = length_label_letter), 
            fontface = "bold", size = 4, 
            show.legend = FALSE)+
  geom_errorbar(aes(ymin=mean_len-se_len, ymax=mean_len+se_len, group=Treatment),
                width=0.05)+
  theme_bw()+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.background = element_blank())+
  ylab("Length (cm)") +
  xlab("Time (mph)")

line_length

## both 
both<-ggpubr::ggarrange(line_length, line_weight, labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
both

## save plot
ggsave("out/figure_line_both.pdf", both, width=210*0.9, height=297*0.33, units="mm", device="pdf")


