# Returns a table of various Bayesian false discover rates given a vector of marginal
# posterior probabilities of inclusion
bfdr_wrap <- function( mppi_vector ){
  BFDR <- numeric()
  for( i in seq(0.05,0.5,0.05)){
    okay <- bfdr(mppi_vector,i)
    BFDR <- rbind( BFDR, c(i,okay$threshold,sum(okay$selected)))
    colnames( BFDR ) <- c( "FDR", "Threshold", "Selected")
  }
  return( BFDR )
}