\name{clean_pathways}
\alias{clean_pathways}
\title{
Clean Enriched Pathways
}
\description{
Clean set of enriched pathways obtained from 'ConsensusPathDB' for gene set overexpression analysis. We also append two columns indicating the number of up-regulated and number of down-regulated genes based on fold change information available in data frame case_vs_ctrl. We cluster pathways based on similarity of gene set using igraph's walktrap clustering algorithm. Within each cluster, pathways are ordered by most to least significant pathway in terms of p-values.
}
\usage{
clean_pathways(df_case_vs_ctrl, df_pathway)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_case_vs_ctrl}{
  Data frame which has at least 2 columns: <gene,fc>.
  Here gene represents the set of genes which are differentially expressed between case and control.
  Here fc represents the fold-change value for each gene.}
  \item{df_pathway}{
    The tab-separated data frame with the pathways information obtained after performing gene set overexpression analysis using 'ConsensusPathDB'.
  }
}
\value{
Returns clean enriched pathways data frame. The data frame has an additional column clusters highlighting the cluster to which each enriched pathway belongs.
}
\author{
Raghvendra Mall
}
\note{
rmall@hbku.edu.qa
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{clean_go_terms}}, \code{\link{clean_pc}}
}
\examples{
data("t.tests.treatment.sign")
data("enriched_pathways")
revised_pathway <- clean_pathways(df_case_vs_ctrl=t.tests.treatment.sign,
                                  df_pathway = enriched_pathways)
print(head(revised_pathway))

}