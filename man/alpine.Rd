\name{alpine}
\docType{data}
\alias{alpine}
\title{Alpine plant community data from the Swiss Alps}
\description{
Presence-absence records of 175 vascular plant species at 912 sites across an elevational gradient (410--3060 m) in the Swiss Alps, with eight environmental predictors.
}
\usage{
data(alpine)
}
\format{
A list with two components:
\describe{
  \item{Y}{A binary matrix of 912 sites by 175 species. Entries are 0 (absent) or 1 (present).}
  \item{X}{A data frame of 912 sites by 8 environmental predictors:
    \itemize{
      \item \code{X}: Swiss coordinate (east), metres.
      \item \code{Y}: Swiss coordinate (north), metres.
      \item \code{DDEG0}: Degree days above 0\eqn{^\circ}C.
      \item \code{SLOPE}: Terrain slope, degrees.
      \item \code{MIND}: Moisture index.
      \item \code{SOLRAD}: Annual solar radiation.
      \item \code{TPI}: Topographic position index.
      \item \code{ELEV}: Elevation, metres (range 410--3060).
    }
  }
}
}
\details{
Presence-absence data for 175 vascular plant species sampled at 912 relevé plots distributed across a broad elevational gradient in the Swiss Alps. Species composition and occurrence probabilities vary strongly with elevation and other environmental predictors, making this dataset suitable for demonstrating joint species distribution models, ordination, species richness prediction, and trait-based analyses.

The species data were originally compiled by D'Amen et al. (2017) and are archived on Dryad (\doi{10.5061/dryad.8mv11}). Elevation values were updated by van der Veen et al. (2021) and are available at Dryad (\doi{10.5061/dryad.pnvx0k6m1}).
}
\references{
D'Amen, M., Mod, H. K., Gotelli, N. J. and Guisan, A. (2018). Disentangling biotic interactions, environmental filters, and dispersal limitation as drivers of species co-occurrence. \emph{Ecography}, 41(8), 1233--1244. \doi{10.1111/ecog.03148}

van der Veen, B., Hui, F. K. C., Hovstad, K. A., Solbu, E. B. and O'Hara, R. B. (2021). Model-based ordination for species with unequal niche widths. \emph{Methods in Ecology and Evolution}, 12(7), 1288--1300. \doi{10.1111/2041-210X.13595}
}
\examples{
data(alpine)
Y <- alpine$Y
X <- alpine$X

\dontrun{
# Fit a JSDM with elevation and two latent variables
X_scaled <- as.data.frame(scale(X[, c("ELEV", "SLOPE", "DDEG0")]))
model <- gllvm(y = Y, X = X_scaled,
               formula = ~ ELEV + I(ELEV^2),
               num.lv = 2,
               family = "binomial")
}
}
