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

clean_go_terms <- function(df_case_vs_ctrl, df_goterms)
{
  #Gene gene names as character and FC into binary based on >1 cutoff
  df_case_vs_ctrl$Genes <- as.character(as.vector(df_case_vs_ctrl$gene))
  df_case_vs_ctrl$Status <- as.numeric(as.vector(df_case_vs_ctrl$fc>1))

  #Convert factors into strings
  df_goterms$term_goid <- as.character(as.vector(df_goterms$term_goid));
  df_goterms$term_category <- as.character(as.vector(df_goterms$term_category));
  df_goterms$term_name <- as.character(as.vector(df_goterms$term_name))
  
  #Find unique GOTERMs
  indices <- which(df_goterms$term_goid %in% unique(df_goterms$term_goid));
  revised_df_goterms <- df_goterms[indices,];
  revised_df_goterms$members_input_overlap_geneids <- as.character(as.vector(revised_df_goterms$members_input_overlap_geneids))
  goterm_genelists <- as.character(as.vector(revised_df_goterms$members_input_overlap_geneids));
  
  #Add the no of up and down regulated genes
  genes_up <- NULL;
  genes_down <- NULL
  for (i in 1:length(goterm_genelists))
  {
    genes_in_goterm <- unlist(strsplit(goterm_genelists[i],split="; "));
    var_1 <- sum(unique(df_case_vs_ctrl[df_case_vs_ctrl$Genes %in% genes_in_goterm,])$Status==1);
    var_2 <- sum(unique(df_case_vs_ctrl[df_case_vs_ctrl$Genes %in% genes_in_goterm,])$Status==0);
    genes_up <- c(genes_up,var_1);
    genes_down <- c(genes_down,var_2);
  }
  revised_df_goterms$genes_up <- genes_up;
  revised_df_goterms$genes_down <- genes_down;

  return(revised_df_goterms)
}
                           
