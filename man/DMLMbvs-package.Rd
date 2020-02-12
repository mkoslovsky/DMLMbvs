\name{DMLMbvs-package}
\alias{DMLMbvs-package}
\alias{DMLMbvs}
\docType{package}
\title{
\packageTitle{DMLMbvs}
}
\description{
Performs Bayesian variable selection for multivariate count data to predict continuous response while simultaneously identifying factors associated with multivariate count data 
}
\usage{
dm_lm_bvs_R()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iterations}{Number of MCMC samples}
  \item{thin}{Thin MCMC chain by # thin}
  \item{y}{subject x 1 vector of continuous response}
  \item{z}{subject x part matrix of multivariate count data}
  \item{x}{subject x covariate matrix of measures}
  \item{alpha}{part X 1 vector of initial intercept values}
  \item{phi}{part X covariate matrix of initial regression coefficients}
  \item{psi}{subject X part matrix of initial compositional probabilities}
  \item{zeta}{part X covariate matrix of initial inclusion indicators for covariates ( 0 or 1 )}
  \item{xi}{part - 1 X 1 vector of initial inclusion indicators for balances}
  \item{sigma2_alpha}{prior value for alpha variance }
  \item{sigma2_phi}{prior value for phi variance}
  \item{h_alpha}{prior variance for intercept term a_0}
  \item{h_beta}{prior variance for beta terms for balances}
  \item{a_m}{parameter for beta prior for balance inclusion probability}
  \item{b_m}{parameter for beta prior for balance inclusion probability}
  \item{a}{parameter for beta prior for covariate inclusion probability}
  \item{b}{parameter for beta prior for covariate inclusion probability}
  \item{a_0}{shape parameter for inverse-gamma prior for sigma}
  \item{b_0}{scale parameter for inverse-gamme prior for sigma }
  \item{seed}{set random seed for simulated data }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Koslovsky, M.D., Hoffman, K., Daniel-MacDougall, C., Vannucci, M.  A Bayesian model of microbiome data for simultaneous identification of covariate associations and prediction of phenotypic outcomes. 2018+.
}
\author{
Matthew Koslovsky email: mkoslovsky12@gmail.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data <- simulate_DMLM()
run <- dm_lm_bvs_R( y = data$Y, z = data$Z, x = data$X )
output <- selected( run, plotting = TRUE )

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bayesian variable selection}
\keyword{Multivariate count data}
\keyword{Compositional Data}
\keyword{Stochastic search variable selection}
\keyword{Microbiome}