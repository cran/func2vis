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

plot_pathways <- function(final_df_pathway, total_no_background_genes=22400, fontsize=9)
{
  #Manipulate Final Pathway data frame
  output_df_pathway <- final_df_pathway[,c("p.value","q.value","pathway","external_id","members_input_overlap","members_input_overlap_geneids","size","effective_size","clusters")]
  colnames(output_df_pathway) <- c("pvalue","qvalue","Description","ID","genes","geneID","size","effective_size","clusters")
  GeneRatio <- NULL
  generatio <- NULL
  BgRatio <- NULL
  bgratio <- NULL
  Count <- NULL
  for (i in 1:nrow(output_df_pathway))
  {
    new_geneID <- paste0(unlist(strsplit(output_df_pathway[i,]$geneID,split="; ")),collapse="/")
    Count <- c(Count,length(unlist(strsplit(output_df_pathway[i,]$geneID,split="; "))))
    GeneRatio <- c(GeneRatio,paste0(length(unlist(strsplit(output_df_pathway[i,]$geneID,split="; "))),"/",output_df_pathway[i,]$size))
    BgRatio <- c(BgRatio,paste0(output_df_pathway[i,]$effective_size,"/",total_no_background_genes))
    generatio <- c(generatio, length(unlist(strsplit(output_df_pathway[i,]$geneID,split="; ")))/output_df_pathway[i,]$size)
    bgratio <- c(bgratio,output_df_pathway[i,]$effective_size/total_no_background_genes)
    output_df_pathway[i,]$geneID <- new_geneID
  }
  output_df_pathway$Count <- as.numeric(as.vector(Count))
  output_df_pathway$GeneRatio <- as.character(as.vector(GeneRatio))
  output_df_pathway$BgRatio <- as.character(as.vector(BgRatio))
  output_df_pathway$generatio <- as.numeric(as.vector(generatio))
  output_df_pathway$bgratio <- as.numeric(as.vector(bgratio))
  output_df_pathway <- output_df_pathway[,c(4,3,11,12,1,2,5,6,10,9,13,14)]
  
  #Color and order the pathways based on generatio
  #colors <- c('#3cb44b', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#46f0f0' , '#bcf60c', '#a9a9a9', '#e6194b', '#ffe119', '#ffe119', '#e6beff')
  #color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  #colors = sample(color, length(unique(output_df_pathway$clusters)))
  n <- length(unique(output_df_pathway$clusters))
  colors <- distinctColorPalette(n)
  output_df_pathway$Description <- paste0(output_df_pathway$Description," [ ",output_df_pathway$GeneRatio," ] ")
  unique_generatios <- unique(output_df_pathway$generatio)
  
  #Make generatios unique
  for (i in 1:length(unique_generatios))
  {
    ids <- which(output_df_pathway$generatio==unique_generatios[i])
    for (j in 1:length(ids))
    {
      output_df_pathway[ids[j],]$generatio <- output_df_pathway[ids[j],]$generatio+j*1e-5
    }
  }
  output_df_pathway <- output_df_pathway[order(output_df_pathway$generatio),]
  output_df_pathway$clusters <- as.factor(output_df_pathway$clusters)
  generatio <- output_df_pathway$generatio
  Description <- output_df_pathway$Description
  pvalue <- output_df_pathway$pvalue
    
  #Plot the pathways
  p <- ggplot(data=output_df_pathway,aes(x=generatio,
                                          y=reorder(Description,generatio),
                                          size=-log10(pvalue)))+
    geom_point(aes(color=-log10(pvalue)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
    scale_size_continuous(name = "-log10(P.adjust)")+   guides(color=guide_legend(), size = guide_legend())+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle("Enriched Pathways") + theme(text = element_text(size=16,color="black")) + 
    theme(axis.text = element_text(size=fontsize, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) +
    theme(axis.text.y=element_text(colour=colors[output_df_pathway$clusters]))
  
  return(p)
}