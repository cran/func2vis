\name{clean_go_terms}
\alias{clean_go_terms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clean Gene Ontologies (GO) Terms
}
\description{
Clean set of enriched goterms obtained from 'ConsensusPathDB' for gene set overexpression analysis. We also append two columns indicating the number of up-regulated and number of down-regulated genes based on fold change information available in data frame case_vs_ctrl.
}
\usage{
clean_go_terms(df_case_vs_ctrl, df_goterms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_case_vs_ctrl}{
  Data frame which has at least 2 columns: <gene,fc>.
  Here gene represents the set of genes which are differentially expressed between case and control.
  Here fc represents the fold-change value for each gene.
}
  \item{df_goterms}{
  The tab-separated data frame with the goterms information obtained after performing gene set overexpression analysis using 'ConsensusPathDB'.
}
}
\value{
Returns clean enriched GO terms data frame.
}
\author{
Raghvendra Mall
}
\note{
rmall@hbku.edu.qa
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See Also as \code{\link{clean_pc}}, \code{\link{plot_go_terms}}
}
\examples{
data("t.tests.treatment.sign")
data("enriched_goterms")
revised_goterms <- clean_go_terms(df_case_vs_ctrl=t.tests.treatment.sign,
                                  df_goterms = enriched_goterms)
print(head(revised_goterms))
}