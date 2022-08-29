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

clean_pc <- function(df_case_vs_ctrl,df_pc)
{
  
  #Gene gene names as character and FC into binary based on >1 cutoff
  df_case_vs_ctrl$Genes <- as.character(as.vector(df_case_vs_ctrl$gene))
  df_case_vs_ctrl$Status <- as.numeric(as.vector(df_case_vs_ctrl$fc>0))
  
  #Convert factors to string variables
  df_pc$complex_name <- as.character(as.vector(df_pc$complex_name));
  df_pc$source <- as.character(as.vector(df_pc$source));
  df_pc$members_input_overlap <- as.character(as.vector(df_pc$members_input_overlap))
  
  # Find unique Protein Complexes by matching all the gene list and removing similar PCs (only keep one PC)
  remove_id <- NULL;
  genes_up <- NULL;
  genes_down <- NULL;
  pc_genelists <- df_pc$members_input_overlap
  for (i in 1:length(pc_genelists))
  {
    genes_in_pc_i <- unlist(strsplit(pc_genelists[i],split="; "));
    var_1 <- sum(unique(df_case_vs_ctrl[df_case_vs_ctrl$Genes %in% genes_in_pc_i,])$Status==1);
    var_2 <- sum(unique(df_case_vs_ctrl[df_case_vs_ctrl$Genes %in% genes_in_pc_i,])$Status==0);
    genes_up <- c(genes_up,var_1);
    genes_down <- c(genes_down,var_2);
    #Match 2 PC and if have similar gene set then keep only 1
    for (j in i:length(pc_genelists))
    {
      if (i!=j)
      {
        genes_in_pc_j <- unlist(strsplit(pc_genelists[j],split="; "))
        if (length(genes_in_pc_i)==length(genes_in_pc_j))
        {
          if (sum(genes_in_pc_i %in% genes_in_pc_j)==length(genes_in_pc_i))
          {
            remove_id <- c(remove_id,j);
          }
        }
      }
    }
  }
  df_pc$genes_up <- genes_up;
  df_pc$genes_down <- genes_down;
  
  if (!is.null(remove_id)){
    revised_df_pc <- df_pc[-remove_id,];
  } else{
    revised_df_pc <- df_pc;
  }
  return(revised_df_pc)
}