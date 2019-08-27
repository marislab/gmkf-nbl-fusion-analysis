# Plot of unique fusion calls per GeneA/B per chromosome per variable
# Only plot in-frame and frame-shift fusions 

library(ggplot2)
source('~/Projects/DGD_Mendelian_RNASeq/R/pubTheme.R')

chrom.plot <- function(fusion.calls, var = "risk", filepath){
  
  fusion.calls <- fusion.calls[grep("other", fusion.calls$Fusion_Type, invert = TRUE),]
  fusion.calls$GeneA_chr <- gsub("[:].*", "", fusion.calls$GeneA_bp)
  fusion.calls$GeneA_chr <- factor(fusion.calls$GeneA_chr, levels = c(1:22,"X","Y"))
  fusion.calls$GeneB_chr <- gsub("[:].*", "", fusion.calls$GeneB_bp)
  fusion.calls$GeneB_chr <- factor(fusion.calls$GeneB_chr, levels = c(1:22,"X","Y"))
  
  # format some naming conventions
  fusion.calls$Caller[fusion.calls$Caller == "arriba"] <- "Arriba"
  fusion.calls$risk <- stringr::str_to_title(fusion.calls$risk)
  
  # split by event type
  inter.chrom <- fusion.calls[which(fusion.calls$GeneA_chr != fusion.calls$GeneB_chr),]
  intra.chrom <- fusion.calls[which(fusion.calls$GeneA_chr == fusion.calls$GeneB_chr),]
  
  # plot intra chrom events
  intra.ct <- plyr::count(intra.chrom, c('Caller','risk','GeneA_chr'))
  
  p <- ggplot(intra.ct, aes(x = GeneA_chr, y = freq, fill = risk)) + 
    geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
    facet_wrap(~Caller, nrow = 2, scales = "free_y") +
    theme_Publication2() + ggtitle("Intra-chromosomal events") + xlab("Chromosome") + 
    labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal")
  
  # plot inter chrom events
  inter.a.ct <- plyr::count(inter.chrom, c('Caller','risk','GeneA_chr'))
  inter.a.ct$Position <- "5'"
  colnames(inter.a.ct)[3] <- "Chr"
  inter.b.ct <- plyr::count(inter.chrom, c('Caller','risk','GeneB_chr'))
  inter.b.ct$Position <- "3'"
  colnames(inter.b.ct)[3] <- "Chr"
  inter.ct <- rbind(inter.a.ct, inter.b.ct)
  inter.ct$Position <- factor(inter.ct$Position, levels = c("5'","3'"))
  
  q <- ggplot(inter.ct, aes(x = Chr, y = freq, fill = risk)) + 
    geom_bar(stat = 'identity', position = 'dodge', color = "black", size = 0.2) +
    facet_grid(Caller~Position, scales = "free", switch = "y") +
    theme_Publication2() + ggtitle("Inter-chromosomal events") + xlab("Chromosome") + 
    labs(fill = "Risk") + theme(legend.position = 'bottom', legend.direction = "horizontal")
  
  # save plots
  ggsave(filename = paste0(filepath, '/Intra-chrom-events.pdf'), plot = p, device = "pdf", width = 8, height = 6)
  ggsave(filename = paste0(filepath, '/Inter-chrom-events.pdf'), plot = q, device = "pdf", width = 12, height = 6)
}
