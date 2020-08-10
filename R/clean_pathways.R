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
clean_pathways <- function(df_case_vs_ctrl, df_pathway)
{
  #Gene gene names as character and FC into binary based on >1 cutoff
  df_case_vs_ctrl$Genes <- as.character(as.vector(df_case_vs_ctrl$gene))
  df_case_vs_ctrl$Status <- as.numeric(as.vector(df_case_vs_ctrl$fc>1))
  
  # Convert factors to string 
  df_pathway$pathway <- as.character(as.vector(df_pathway$pathway));
  df_pathway$source <- as.character(as.vector(df_pathway$source));
  df_pathway$members_input_overlap <- as.character(as.vector(df_pathway$members_input_overlap))
  df_pathway$external_id <- as.character(as.vector(df_pathway$external_id))
  df_pathway$members_input_overlap_geneids <- as.character(as.vector(df_pathway$members_input_overlap_geneids))
  
  # Find unique PATHWAYs by matching all the gene list
  remove_id <- NULL;
  genes_up <- NULL;
  genes_down <- NULL;
  pathway_genelists <- df_pathway$members_input_overlap
  for (i in 1:length(pathway_genelists))
  {
    genes_in_pathway_i <- unlist(strsplit(pathway_genelists[i],split="; "));
    var_1 <- sum(unique(df_case_vs_ctrl[df_case_vs_ctrl$Genes %in% genes_in_pathway_i,])$Status==1);
    var_2 <- sum(unique(df_case_vs_ctrl[df_case_vs_ctrl$Genes %in% genes_in_pathway_i,])$Status==0);
    genes_up <- c(genes_up,var_1);
    genes_down <- c(genes_down,var_2);
    #Match 2 Pathways and if have similar gene set then keep only 1
    for (j in i:length(pathway_genelists))
    {
      if (i!=j)
      {
        genes_in_pathway_j <- unlist(strsplit(pathway_genelists[j],split="; "))
        #Two pathways with exact same gene list, keep one
        if (length(genes_in_pathway_i)<=length(genes_in_pathway_j))
        {
          if (sum(genes_in_pathway_i %in% genes_in_pathway_j)==length(genes_in_pathway_i))
          {
            remove_id <- c(remove_id,j);
          }
        }
        if (df_pathway$pathway[i]==df_pathway$pathway[j])
        {
          remove_id <- c(remove_id,j);
        }
      }
    }
  }
  df_pathway$genes_up <- genes_up;
  df_pathway$genes_down <- genes_down;
  
  if (!is.null(remove_id)){
    revised_df_pathway <- df_pathway[-unique(remove_id),];
  } else{
    revised_df_pathway <- df_pathway;
  }
  
  #Cluster the pathways based on overlap in gene lists in pathways
  pathway_genelists <- revised_df_pathway$members_input_overlap;
  pathway_pathway_similarity_matrix <- matrix(0,nrow=length(pathway_genelists),ncol=length(pathway_genelists));
  for (i in 1:length(pathway_genelists))
  {
    for (j in 1:length(pathway_genelists))
    {
      if (i==j)
      {
        pathway_pathway_similarity_matrix[i,j] <- 1;
      } else {
        genes_in_pathway_i <- unlist(strsplit(pathway_genelists[i],split="; "))
        genes_in_pathway_j <- unlist(strsplit(pathway_genelists[j],split="; "))
        common_genes <- intersect(genes_in_pathway_i,genes_in_pathway_j);
        all_genes <- union(genes_in_pathway_i,genes_in_pathway_j);
        pathway_pathway_similarity_matrix[i,j] <- length(common_genes)/length(all_genes)
      }
    }
  }
  g <- graph_from_adjacency_matrix(pathway_pathway_similarity_matrix,mode=c("directed"),weighted = TRUE)
  comm_model <- walktrap.community(g,steps=4,merges = TRUE)
  comm_memb <- membership(comm_model)
  revised_df_pathway$clusters <- comm_memb;
  
  #Order pathways in all clusters in increasing order of pvalue
  revised_df_pathway <- revised_df_pathway[order(revised_df_pathway$clusters,decreasing = F),];
  unique_clusters <- unique(revised_df_pathway$clusters)
  
  #Order all pathways in a cluster based on significance
  for (i in 1:length(unique_clusters))
  {
    clusterid <- unique_clusters[i];
    indices <- which(revised_df_pathway$clusters==clusterid)
    if (length(indices)>1)
    {
      revised_df_pathway[indices,]  <- revised_df_pathway[indices[order(revised_df_pathway[indices,]$p.value,decreasing=F)],]
    }
  }
  
  #Order all pathways based on increasing pvalue
  final_df_pathway <- NULL;
  old_revised_df_pathway <- revised_df_pathway;
  while(nrow(revised_df_pathway)>0)
  {
    min_pvalue <- min(revised_df_pathway$p.value);
    clusterid <- revised_df_pathway[revised_df_pathway$p.value==min_pvalue,]$clusters;
    if (length(clusterid)==1)
    {
      ids <- which(revised_df_pathway$clusters==clusterid);
      final_df_pathway <- rbind(final_df_pathway,old_revised_df_pathway[old_revised_df_pathway$clusters==clusterid,]);
      revised_df_pathway <- revised_df_pathway[-ids,];
    } else {
      for (j in 1:length(clusterid))
      {
        ids <- which(revised_df_pathway$clusters==clusterid[j]);
        final_df_pathway <- rbind(final_df_pathway,old_revised_df_pathway[old_revised_df_pathway$clusters==clusterid[j],]);
        revised_df_pathway <- revised_df_pathway[-ids,];
      }
    }
  }
  
  #Get characteristics of all the columns in final pathway dataframe
  final_df_pathway <- as.data.frame(final_df_pathway)
  final_df_pathway$p.value <- as.numeric(as.vector(final_df_pathway$p.value))
  final_df_pathway$q.value <- as.numeric(as.vector(final_df_pathway$q.value))
  final_df_pathway$pathway <- as.character(as.vector(final_df_pathway$pathway))
  final_df_pathway$source <- as.character(as.vector(final_df_pathway$source))
  final_df_pathway$external_id <- as.character(as.vector(final_df_pathway$external_id))
  final_df_pathway$members_input_overlap <- as.character(as.vector(final_df_pathway$members_input_overlap))
  final_df_pathway$members_input_overlap_geneids <- as.character(as.vector(final_df_pathway$members_input_overlap_geneids))
  final_df_pathway$size <- as.numeric(as.vector(final_df_pathway$size))
  final_df_pathway$effective_size <- as.numeric(as.vector(final_df_pathway$effective_size))
  final_df_pathway$genes_up <- as.numeric(as.vector(final_df_pathway$genes_up))
  final_df_pathway$genes_down <- as.numeric(as.vector(final_df_pathway$genes_down))
  final_df_pathway$clusters  <- as.numeric(as.vector(final_df_pathway$clusters))

  return(final_df_pathway)
}