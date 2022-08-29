\name{plot_go_terms}
\alias{plot_go_terms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bupple Plot for GO Terms
}
\description{
Make a bubble plot for significantly enriched Gene Ontologies (GO) Terms obtained after performing gene set overexpression analysis using 'ConsensusPathDB'. 
}
\usage{
plot_go_terms(df_goterms, total_no_background_genes, negative_log_10_p_value_cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_goterms}{
   The tab-separated data frame with the GO terms information obtained after performing gene set overexpression analysis   using 'ConsensusPathDB'.
   }
  \item{total_no_background_genes}{
  Total no of genes in the background set.
  }
  \item{negative_log_10_p_value_cutoff}{
  The threshold on -log10(pvalue) to be used to identify the GO terms to be highlighted in the plot.
}
}
\details{
Plots the significantly enriched molecular function (m), cellular components (c) and biological processes (b) obtained via ConsensusPathDB.
}
\value{
Returns a bubble plot of type ggplot.
}
\author{
Raghvendra Mall
}
\note{
rmall@hbku.edu.qa
}

\examples{
data("enriched_goterms")
g <- plot_go_terms(df_goterms = enriched_goterms, negative_log_10_p_value_cutoff=17)
g
}