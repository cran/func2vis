\name{clean_pc}
\alias{clean_pc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clean Enriched Protein Complexes
}
\description{
Clean set of enriched protein complexes obtained from 'ConsensusPathDB' for gene set overexpression analysis. We also append two columns indicating the number of up-regulated and number of down-regulated genes based on fold change information available in data frame case_vs_ctrl.
}
\usage{
clean_pc(df_case_vs_ctrl,df_pc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_case_vs_ctrl}{
  Data frame which has at least 2 columns: <gene,fc>.
  Here gene represents the set of genes which are differentially expressed between case and control.
  Here fc represents the fold-change value for each gene.
  }
  \item{df_pc}{
     The tab-separated data frame with the protein complexes information obtained after performing gene set overexpression analysis using 'ConsensusPathDB'.
  }
}

\value{
Returns clean enriched protein complexes data frame.
}
\author{
Raghvendra Mall
}
\note{
rmall@hbku.edu.qa
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{clean_go_terms}}, \code{\link{plot_go_terms}}
}
\examples{
data("t.tests.treatment.sign")
data("enriched_pc")
revised_pc <- clean_pc(df_case_vs_ctrl=t.tests.treatment.sign,
                                  df_pc = enriched_pc)
print(head(revised_pc))
}

