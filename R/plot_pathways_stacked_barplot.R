#    This is an implementation of func2vis algorithm for visualization and
#    cleanup of enriched GO Terms, Protein Complexes and Pathways obtained
#    using gene set overexpression analysis using ConsensusPathDB.
#    Copyright (C) 2020  Raghvendra Mall

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program, see LICENSE.
plot_pathways_stacked_barplot <- function(final_df_pathway)
{
  #Make dataframe for stacked barplot
  plot_df <- final_df_pathway
  plot_df$p.value <- round(plot_df$p.value,7)
  for (i in 1:nrow(plot_df))
  {
    plot_df$pathway[i] <- paste0(plot_df$pathway[i]," x [ ",plot_df$genes_up[i]+plot_df$genes_down[i],"/",plot_df$effective_size[i]," ]");
  }
  final_plot_df <- NULL
  #Get percentage of up, down and non-differentially expressed genes in a pathway
  for (i in 1:nrow(plot_df))
  {
    var1 <- cbind(plot_df$pathway[i],round((plot_df$genes_up[i]/plot_df$effective_size[i])*100,0),plot_df$genes_up[i],plot_df$p.value[i],"Genes Up",plot_df$clusters[i])
    var2 <- cbind(plot_df$pathway[i],round((plot_df$genes_down[i]/plot_df$effective_size[i])*100,0),plot_df$genes_down[i],plot_df$p.value[i],"Genes Down",plot_df$clusters[i])
    var3 <- cbind(plot_df$pathway[i],100-round((plot_df$genes_up[i]/plot_df$effective_size[i])*100,0)-round((plot_df$genes_down[i]/plot_df$effective_size[i])*100,0),
                  plot_df$effective_size[i],plot_df$p.value[i],"All Other Genes",plot_df$clusters[i])
    final_plot_df <- rbind(final_plot_df,var1);
    final_plot_df <- rbind(final_plot_df,var2);
    final_plot_df <- rbind(final_plot_df,var3);
  }
  final_plot_df <- as.data.frame(final_plot_df)
  colnames(final_plot_df) <- c("Pathway","Genes_Type","Genes_Number","Pvalue","Gene_Group","Cluster")
  final_plot_df$Pathway <- as.character(as.vector(final_plot_df$Pathway))
  final_plot_df$Genes_Type <- as.numeric(as.vector(final_plot_df$Genes_Type))
  final_plot_df$Genes_Number <- as.numeric(as.vector(final_plot_df$Genes_Number))
  final_plot_df$Pvalue <- as.character(as.vector(final_plot_df$Pvalue))
  final_plot_df$Gene_Group <- as.character(as.vector(final_plot_df$Gene_Group))
  final_plot_df$Cluster <- as.numeric(as.vector(final_plot_df$Cluster))
  final_plot_df$Pathways <- factor(final_plot_df$Pathway,levels=unique(plot_df$pathway))
  pathwaylabs <- plot_df$pathway;
  pvaluelabs <- factor(round(-log10(plot_df$p.value+.Machine$double.eps),2),levels=unique(round(-log10(plot_df$p.value+.Machine$double.eps),2)))
  cluster <- plot_df$clusters;
  n <- length(unique(cluster))
  color_palette = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  color_palette = color_palette[-grep("white",color_palette)]
  color_palette = color_palette[sample(length(color_palette))]
  palette <- color_palette[1:n]
  Pathways <- final_plot_df$Pathways
  Genes_Type <- final_plot_df$Genes_Type
  Gene_Group <- final_plot_df$Gene_Group
  
  #Make stacked barplot
  p <- ggplot(data=final_plot_df, aes(x=rev(as.numeric(Pathways)),y=Genes_Type,fill=Gene_Group,label=Genes_Type)) +
    geom_col(position = position_stack()) + coord_flip() +
    scale_fill_manual(values=c("Grey","Green","Red")) +
    theme_minimal() + ylab("Percentage of Genes") + xlab("Pathways") +
    scale_x_continuous(breaks=1:length(pathwaylabs),
                       labels=rev(pathwaylabs),
                       sec.axis=sec_axis(~.,name="-log10(Pvalue)",
                                         breaks=1:length(pvaluelabs),
                                         labels=rev(pvaluelabs)),
                       expand = c(0,0.6)) +
    theme(axis.text.y = element_text(hjust = 1, colour = rev(palette[cluster]))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(fill="") + ggtitle("Enriched Pathways for Treatment samples at DE genes p-value = 0.05") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}