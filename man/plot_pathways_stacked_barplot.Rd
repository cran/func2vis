\name{plot_pathways_stacked_barplot}
\alias{plot_pathways_stacked_barplot}
\title{
Stacked Barplot of Cleaned Pathways
}
\description{
Make a stacked barplot like the one available in Ingenuity Pathway Analysis highlighting percentage of up, down and non-differentially expressed genes in the set of clean enriched pathways obtained from 'ConsensusPathDB' by performing gene set overexpression analysis. You need to run the function \code{\link{clean_pathways}} to obtain the input data frame 
}
\usage{
plot_pathways_stacked_barplot(final_df_pathway)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{final_df_pathway}{
  Clean and clustered pathways obtained using \code{\link{clean_pathways}}.
}
}
\value{
Returns a stacked barplot of type ggplot.
}
\author{
Raghvendra Mall
}
\note{
rmall@hbku.edu.qa
}
\seealso{
\code{\link{clean_pathways}}, \code{\link{plot_go_terms}}
}
\examples{
data("t.tests.treatment.sign")
data("enriched_pathways")
revised_pathway <- clean_pathways(df_case_vs_ctrl=t.tests.treatment.sign,
                                  df_pathway = enriched_pathways)
p <- plot_pathways_stacked_barplot(revised_pathway)
p
}