setwd('~/Projects/Maris-lab/GMKF_NBL/')

source('R/piedonut_utils.R')
source('~/Projects/DGD_Mendelian_RNASeq/R/pubTheme.R')

library(ggplot2)
library(moonBook)
library(webr)
library(ggforce)
library(grid)

data("acs")
acs <- acs[,c('Dx','smoking')]

load('data/arriba_merged.RData')
arriba.merge$Fusion_Caller <- "Arriba"
arriba.merge <- arriba.merge[,c("Fusion_Caller","confidence")]
load('data/star_merged.RData')
star.merge$Fusion_Caller <- "STAR-Fusion"
star.merge$confidence <- "high"
star.merge <- star.merge[,c("Fusion_Caller","confidence")]
total <- rbind(star.merge, arriba.merge)

total$confidence <- factor(total$confidence, levels = c("low","medium","high"))
PieDonut(total, aes(pies = Fusion_Caller, donuts = confidence), showRatioPie = T, color = c("black"))

ct <- plyr::count(total, c("Fusion_Caller","confidence"))
total <- merge(total, ct, by = c("Fusion_Caller","confidence"))
total <- as.data.frame(total)
total <- unique(total)

ggplot(total, aes(Fusion_Caller, y = freq, fill = confidence, label = freq)) + 
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("high"= "#D11141",
                               "medium"="#F37735",
                               "low"="#FFC425")) + theme_bw() + theme_Publication2()

