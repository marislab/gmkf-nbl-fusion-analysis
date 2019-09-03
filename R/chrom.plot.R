# Plot of unique fusion calls per GeneA/B per chromosome per variable
# Only plot in-frame and frame-shift fusions 

library(ggplot2)
library(GenomicRanges)
library(ggpubr)
source('~/Projects/DGD_Mendelian_RNASeq/R/pubTheme.R')

chrom.plot <- function(fusion.calls, var = "risk", filepath, cytoband){
  
  fusion.calls <- fusion.calls[grep("other", fusion.calls$Fusion_Type, invert = TRUE),]
  fusion.calls$GeneA_chr <- gsub("[:].*", "", fusion.calls$GeneA_bp)
  fusion.calls$GeneA_chr <- factor(fusion.calls$GeneA_chr, levels = c(1:22,"X","Y"))
  fusion.calls$GeneB_chr <- gsub("[:].*", "", fusion.calls$GeneB_bp)
  fusion.calls$GeneB_chr <- factor(fusion.calls$GeneB_chr, levels = c(1:22,"X","Y"))
  
  # merge cytoband info
  fusion.calls <- merge.cytoband.info(fusion.calls, cytoband)
  
  # format some naming conventions
  fusion.calls$Caller[fusion.calls$Caller == "arriba"] <- "Arriba"
  fusion.calls$risk <- stringr::str_to_title(fusion.calls$risk)
  
  # split by event type
  inter.chrom <- fusion.calls[which(fusion.calls$GeneA_chr != fusion.calls$GeneB_chr),]
  intra.chrom <- fusion.calls[which(fusion.calls$GeneA_chr == fusion.calls$GeneB_chr),]
  
  # plot intra chrom events
  tmp <- intra.chrom[,c("Caller","risk","GeneA_chr","GeneA_cyto")]
  tmp <- dcast(tmp, Caller+risk+GeneA_chr~GeneA_cyto, fun.aggregate = length, drop = FALSE)
  intra.ct <- melt(tmp, variable.name = "Cyto", value.name = "freq")
  
  p <- ggplot(intra.ct, aes(x = Cyto, y = freq, fill = risk)) +
    geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
    facet_wrap(~Caller, nrow = 2, scales = "free_y") +
    theme_Publication2() + ggtitle("Intra-chromosomal events\n") + xlab("Chromosome") +
    ylab("Frequency of Fusion Events") +
    labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal") +
    theme(plot.title = element_text(face = "bold", size = 14)) + scale_fill_manual(values = c("High"= "#D11141",
                                                                                              "Intermediate"="#F37735",
                                                                                              "Low"="#FFC425"))
  ggsave(filename = paste0(filepath, '/Intra-chrom-events.pdf'), plot = p, device = "pdf", width = 15, height = 6)
  
  # tmp1 <- intra.chrom[,c("Caller","risk","GeneA_chr","GeneA_cyto")]
  # tmp1$Position <- "5'"
  # tmp2 <- intra.chrom[,c("Caller","risk","GeneB_chr","GeneB_cyto")]
  # tmp2$Position <- "3'"
  # colnames(tmp2)[3:4] <- c("GeneA_chr","GeneA_cyto")
  # tmp <- rbind(tmp1, tmp2)
  # tmp <- dcast(tmp, Caller+risk+GeneA_chr+Position~GeneA_cyto, fun.aggregate = length, drop = FALSE)
  # intra.ct <- melt(tmp, variable.name = "Cyto", value.name = "freq")
  # 
  # fig <- ggpubr::ggarrange(
  #   ggplot(intra.ct[which(intra.ct$Caller == "Arriba"),], aes(x = as.factor(Cyto), y = freq, fill = risk)) + 
  #     geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
  #     facet_grid(Position~Caller, scales = "free_y", switch = "y") +
  #     theme_Publication2() + theme(axis.title = element_blank(), 
  #                                  axis.text.x = element_blank(), 
  #                                  axis.ticks.x = element_blank()) + ylab("") +
  #     labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal") +
  #     theme(plot.margin = unit(c(1,1,0.5,1), "lines")),
  #   ggplot(intra.ct[which(intra.ct$Caller == "STARFusion"),], aes(x = as.factor(Cyto), y = freq, fill = risk)) + 
  #     geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
  #     facet_grid(Position~Caller, scales = "free_y", switch = "y") +
  #     theme_Publication2() + theme(axis.title.y = element_blank()) + xlab("Chromosome") + ylab("") +
  #     labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal") +
  #     theme(plot.margin = unit(c(1,1,1,1), "lines")), 
  #   common.legend = TRUE, legend = "bottom", nrow = 2, align = "v") 
  # annotate_figure(fig, top = text_grob("Intra-chromosomal events", face = "bold", size = 14),
  #                 left = text_grob("\nFrequency of Fusion Events", face = "bold", rot = 90))
  
  # plot inter chrom events
  tmp1 <- inter.chrom[,c("Caller","risk","GeneA_chr","GeneA_cyto")]
  tmp1$Position <- "5'"
  tmp2 <- inter.chrom[,c("Caller","risk","GeneB_chr","GeneB_cyto")]
  tmp2$Position <- "3'"
  colnames(tmp2)[3:4] <- c("GeneA_chr","GeneA_cyto")
  tmp <- rbind(tmp1, tmp2)
  tmp1 <- dcast(tmp, Caller+risk+GeneA_chr+Position~GeneA_cyto, fun.aggregate = length, drop = FALSE)
  tmp2 <- dcast(tmp, Caller+risk+GeneA_chr~GeneA_cyto, fun.aggregate = length, drop = FALSE)
  inter.ct <- melt(tmp1, variable.name = "Cyto", value.name = "freq")
  inter.ct$Position <- factor(inter.ct$Position, levels = c("5'","3'"))
  inter.ct.2 <- melt(tmp2, variable.name = "Cyto", value.name = "freq")
  
  
  q <- ggplot(inter.ct.2, aes(x = Cyto, y = freq, fill = risk)) +
    geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
    facet_wrap(~Caller, nrow = 2, scales = "free_y") +
    theme_Publication2() + ggtitle("Inter-chromosomal events\n") + xlab("Chromosome") +
    ylab("Frequency of Fusion Events") +
    labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal") +
    theme(plot.title = element_text(face = "bold", size = 14)) + scale_fill_manual(values = c("High"= "#D11141",
                                                                                              "Intermediate"="#F37735",
                                                                                              "Low"="#FFC425"))
  
  ggsave(filename = paste0(filepath, '/Inter-chrom-events.pdf'), plot = q, device = "pdf", width = 15, height = 6)
  
  fig <- ggpubr::ggarrange(
    ggplot(inter.ct[which(inter.ct$Caller == "Arriba"),], aes(x = as.factor(Cyto), y = freq, fill = risk)) + 
      geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
      facet_grid(Position~Caller, scales = "free_y", switch = "y") +
      theme_Publication2() + theme(axis.title = element_blank(), 
                                   axis.text.x = element_blank(), 
                                   axis.ticks.x = element_blank()) + ylab("") +
      labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal") +
      theme(plot.margin = unit(c(1,1,0.5,1), "lines")) + scale_fill_manual(values = c("High"= "#D11141",
                                                                                      "Intermediate"="#F37735",
                                                                                      "Low"="#FFC425")),
    ggplot(inter.ct[which(inter.ct$Caller == "STARFusion"),], aes(x = as.factor(Cyto), y = freq, fill = risk)) + 
      geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
      facet_grid(Position~Caller, scales = "free_y", switch = "y") +
      theme_Publication2() + theme(axis.title.y = element_blank()) + xlab("Chromosome") + ylab("") +
      labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal") +
      theme(plot.margin = unit(c(0.5,1,1,1), "lines")) + scale_fill_manual(values = c("High"= "#D11141",
                                                                                      "Intermediate"="#F37735",
                                                                                      "Low"="#FFC425")), 
    common.legend = TRUE, legend = "bottom", nrow = 2, align = 'v') 
  annotate_figure(fig, top = text_grob("Inter-chromosomal events", face = "bold", size = 14),
                  left = text_grob("\nFrequency of Fusion Events", face = "bold", rot = 90)) %>%
    ggexport(filename = paste0(filepath, '/Inter-chrom-events_split.pdf'), width = 15, height = 8)
  
  return(fusion.calls)
}
