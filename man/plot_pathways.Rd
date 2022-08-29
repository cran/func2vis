\name{plot_pathways}
\alias{plot_pathways}
\title{
Plot clean enriched pathways as a bubble plot
}
\description{
Make a bubble plot of clean enriched pathways obtained from 'ConsensusPathDB' by performing gene set overexpression analysis. Colours represent the clusters to which each pathway belongs. You need to run the function \code{\link{clean_pathways}} to obtain the input data frame.
}
\usage{
plot_pathways(final_df_pathway, total_no_background_genes, fontsize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{final_df_pathway}{
    Clean and clustered pathways obtained using \link{clean_pathways}.
  }
  \item{total_no_background_genes}{
    Total no of genes in the background set.
  }
  \item{fontsize}{
    Font size of the pathways to be displayed on y-axis.
  }
}
\value{
Returns a bubble plot of type ggplot. Colours represent the clusters to which each pathway belongs. 
}
\author{
Raghvendra Mall
}
\note{
rmall@hbku.edu.qa
}


\seealso{
See Also as \code{\link{clean_pathways}}, \code{\link{plot_pathways_stacked_barplot}}, \code{\link{plot_go_terms}}
}
\examples{
data("t.tests.treatment.sign")
data("enriched_pathways")
revised_pathway <- clean_pathways(df_case_vs_ctrl=t.tests.treatment.sign,
                                  df_pathway = enriched_pathways)
p <- plot_pathways(revised_pathway)
p
}