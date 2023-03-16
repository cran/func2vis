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
plot_go_terms <- function(df_goterms, total_no_background_genes=22400, negative_log_10_p_value_cutoff = 3, max_overlap = 20)
{
  #Convert factors to strings
  df_goterms$term_goid <- as.character(as.vector(df_goterms$term_goid));
  df_goterms$term_category <- as.character(as.vector(df_goterms$term_category));
  df_goterms$term_name <- as.character(as.vector(df_goterms$term_name))
  
  #Find unique GOTERMs
  indices <- which(df_goterms$term_goid %in% unique(df_goterms$term_goid));
  revised_df_goterms <- df_goterms[indices,];
  revised_df_goterms$p.value <- as.numeric(as.vector(revised_df_goterms$p.value))
  revised_df_goterms$term_category <- as.character(as.vector(revised_df_goterms$term_category))
  revised_df_goterms$term_name <- as.character(as.vector(revised_df_goterms$term_name))
  revised_df_goterms$members_input_overlap_geneids <- as.character(as.vector(revised_df_goterms$members_input_overlap_geneids))
  goterm_genelists <- as.character(as.vector(revised_df_goterms$members_input_overlap_geneids));
  
  #Prepare data frame for the bubble plot
  generatio <- NULL
  bgratio <- NULL
  count <- NULL
  term_info <- NULL
  for (i in 1:nrow(revised_df_goterms))
  {
    term_info  <- c(term_info,paste0(as.character(revised_df_goterms[i,]$term_category),as.character(revised_df_goterms[i,]$term_level)))
    geneids <- unlist(strsplit(revised_df_goterms[i,]$members_input_overlap_geneids,split="; "))
    count <- c(count,length(geneids))
    generatio <- c(generatio,length(geneids)/revised_df_goterms[i,]$size)
    bgratio <- c(bgratio,revised_df_goterms[i,]$effective_size/total_no_background_genes)
    revised_df_goterms[i,]$members_input_overlap_geneids <- paste0(geneids,collapse="/");
  }
  revised_df_goterms$count <- count
  revised_df_goterms$generatio <- generatio
  revised_df_goterms$bgratio <- bgratio
  revised_df_goterms$term_info <- term_info

  # Data frame for GOTERMS whose pvalues are below a cutoff to be highlighted in the plot
  temp_df <- revised_df_goterms[-log10(revised_df_goterms$p.value)>negative_log_10_p_value_cutoff,]
  max_x = max(revised_df_goterms$generatio*100);
  min_x = min(revised_df_goterms$generatio*100);
  p.value <- revised_df_goterms$p.value
  term_category <- revised_df_goterms$term_category
  term_name <- temp_df$term_name
  
  #Make the ggplot
  g <- ggplot(revised_df_goterms,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) +
    geom_point(aes(color=as.factor(term_category),fill=as.factor(term_category),size=-log10(p.value))) +
    xlab("Overexpressed Genes to Genes in GO Terms (in Percentage)") + ylab("-log10(Pvalue)") +
    geom_hline(yintercept = 1.301, col="black")+
    geom_hline(yintercept = negative_log_10_p_value_cutoff, col = "blue") +
    geom_vline(xintercept = max(0,min_x+10), col = "black") +
    scale_shape_manual(name = "Category", values = c(15,16,17)) +
    scale_color_manual(name = "Category", values =  c("red","green","blue")) +
    scale_fill_manual(name="Category", values = alpha(c("red","green","blue"),0.5))+
    theme_bw() + xlim(c(0.0,max_x+5))+
    geom_label_repel(data=temp_df, aes(label = term_name, x=100*generatio, y=-log10(p.value)),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     max.overlaps = max_overlap,
                     segment.color = 'grey50') +
    theme(text = element_text(size=20)) + theme(plot.title = element_text(hjust = 0.5))
   return(g)
}