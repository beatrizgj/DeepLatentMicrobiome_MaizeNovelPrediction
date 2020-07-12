suppressPackageStartupMessages({
  library(phyloseq)
  library(biomformat)
  library(ggplot2)
  library(plyr)
  library(dplyr)
})

# Author: Beatriz Garcia-Jimenez
build_physeq_object_otu_tax <- function(fotu,ftax){
  otumat <- read.table(fotu,sep='\t',row.names=1,header=TRUE,check.names=FALSE)
  OTU = otu_table(otumat, taxa_are_rows = TRUE)
  
  taxmat <- read.table(ftax,sep='\t',row.names=1,header=TRUE)
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxmat=as.matrix(taxmat)
  TAX = tax_table(taxmat)
  
  phyObj = phyloseq(OTU, TAX) 
  
  return(phyObj) 
}
###

plot_taxa_agg <- function (phy, tax.rank, temp, rain, age){
  #color according to tax_rank
  if (tax.rank != "Species") {
    #group data according to tax rank (as a parameter)
    phy.rank <- phy %>% 
      tax_glom(taxrank = tax.rank)
  }else{
    phy.rank=phy
  }
  
  phy.rank.df=phy.rank %>%
    psmelt()
 
  #plot
  p <- ggplot(phy.rank.df, aes(x = phy.rank.df$Sample, y = Abundance, fill = phy.rank.df[,tax.rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=tax.rank)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    xlab("Samples")+
    ggtitle("Predicted from novel environment features", subtitle=paste("T=",temp,"F, rain=",rain," inches, plant age=",age," weeks",sep="")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title=element_text(size = 20, face = "bold"),
          plot.subtitle=element_text(size = 16),
          axis.title.y=element_text(size=16))
  # To save in .pdf file
  pdf(paste('barplot_predicted_otus_novel_features',temp,rain,age,tax.rank,'.pdf',sep='_'))  
  plot(p)
  dev.off()
  # To save in .png file
  png(paste('barplot_predicted_otus_novel_features',temp,rain,age,tax.rank,'.png',sep='_'))  
  plot(p)
  dev.off()
  # To show online
  plot(p)
}
