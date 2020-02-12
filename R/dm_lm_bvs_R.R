# Wrapper function for the Rcpp code to initiate with defaults and simulate data if requested
dm_lm_bvs_R <- function( iterations = 30000, thin = 10, y = NULL, z = NULL, x = NULL, alpha = NULL, phi = NULL,
                         psi = NULL, zeta = NULL, xi = NULL, sigma2_alpha = sqrt( 10 ), sigma2_phi = sqrt( 10 ),
                         h_alpha = 1, h_beta = 1, delta = 0, rho = 5*10^6, a_m = 1, b_m = 1, a = 1, b = 1,
                         a_0 = 2, b_0 = 2, subject_sim = 50, B_sim = 50, covariates_sim = 50,
                         covariance = 0.4, seed = 2, num_cov = 10, num_bal = 5){
  library(mvtnorm)
  library(bvsDMLM)
  library(MCMCpack)
  library(Rcpp)
  library(ggplot2)
  library(RcppArmadillo)
  # iterations - Number of MCMC samples, Default = 30000
  # thin - Then MCMC by # thin, Default = 10
  # y - subject x 1 vector of continuous response
  # z - subject x part matrix of multivariate count data
  # x - subject x covariate matrix of measures
  # alpha - part X 1 vector of initial intercept values
  # phi - part X covariate matrix of initial regression coefficients
  # psi - subject X part matrix of initial compositional probabilities
  # zeta - part X covariate matrix of initial inclusion indicators for covariates ( 0 or 1 )
  # xi - part - 1 X 1 vector of initial inclusion indicators for balances
  # sigma2_alpha - prior value for alpha variance, Default = sqrt(10)
  # sigma2_phi - prior value for phi variance, Default = sqrt(10)
  # h_alpha - prior variance for intercept term a_0, Default = 1
  # h_beta - prior variance for beta terms for balances, Default = 1
  # delta - Dirichlet proposal scale parameter, Default = 0
  # rho - Dirichlet proposal shift parameter, Default = 5x10^6
  # a_m -  parameter for beta prior for balance inclusion probability, Default = 1
  # b_m -  parameter for beta prior for balance inclusion probability, Default = 1
  # a -  parameter for beta prior for covariate inclusion probability, Default = 1
  # b -  parameter for beta prior for covariate inclusion probability, Default = 1
  # a_0 - shape parameter for inverse-gamma prior for sigma, Default = 2
  # b_0 - scale parameter for inverse-gamme prior for sigma, Default = 2
  # sumbject_sim - number of subjects to simulate, Default = 50
  # B_sim - number of compositional parts to simulate, Default = 50
  # covariates_sim - number of covariates to simulate, Default = 50
  # covariance - simulated covariance structure, Default = 0.4
  # seed - set random seed for simulated data, Default = 2
  # num_cov - number of associated covariates across compostional parts, Default = 10
  # num_bal - number of assocated balances, Default = 5
  
  # Set seed for replication
  set.seed( seed )
  
  # Simulate data 
  if( is.null( x ) | is.null( y ) | is.null( z ) ){
    
    # Set covariance strucuture for covariates 
    sig <- diag(covariates_sim)
    for( i in 1:covariates_sim ){
      for( j in 1:covariates_sim ){
        if( i != j){
          sig[ i , j] = covariance^abs(i - j)
        }
      }
    }
    
    # Simulate covariates from MVT normal distribution 
    x <- rmvnorm( subject_sim, rep( 0, covariates_sim ), sig )
    
    # Set influential covariates 
    zeta_sim <- matrix( 0, B_sim, covariates_sim)
    true_cov <- cbind( sample( seq( 1,B_sim), num_cov ),sample( seq( 1, covariates_sim ), num_cov) )
    zeta_sim[  true_cov ] <- 1
    
    # Set alpha (intercept) values 
    alpha_sim <- matrix( 1, nrow = subject_sim , ncol = 1)%*%( runif( n = B_sim, -2.3,2.3 ) )
    
    # Set regression coefficients for associated covariates 
    true_coeff <- runif( num_cov, 0.75, 1)
    phi_sim <- matrix( 0, B_sim, covariates_sim)
    phi_sim[ true_cov ] <- true_coeff*sample(c(1,-1), num_cov, replace = TRUE)
    
    # Using covariates and respective regression coefficients, set probabilities for count data 
    inside_sim <- exp( alpha_sim + x%*%t(phi_sim) )
    psi_sim <- inside_sim/rowSums(inside_sim)
    psi_sim_overdispersed = psi_sim*(1-0.01)/0.01
    
    # Simulate probabilities for each branch
    prob_sim <- apply( psi_sim_overdispersed, 1,  function(x){ rdirichlet(1,x) } )
    
    # Simulate count data for each subtree branch (simplified for bifurcating tree in this example)
    z <- t( apply( t( prob_sim ), 1, function(x){ rmultinom(1, sample( seq(2500, 7500), 1 ), x) } ) )
    
    # Adjust count data 
    Z_norm = z
    Z_norm[ Z_norm == 0 ] <- 0.5
    
    Z_norm <- Z_norm/rowSums( Z_norm )
    
    # Calculate Balances
    Balances <- ilr_fun_cpp( Z_norm  )
    
    # Select Balances
    xi_sim <- rep( 0, B_sim - 1 )
    true_cov_beta <- sample( seq( 1,B_sim - 1 ), num_bal ) 
    xi_sim[  true_cov_beta ] <- 1
    
    # Set respective balance coefficients
    true_coeff_beta <- runif( num_bal, 0.75, 1 )
    beta_sim <- rep( 0, B_sim - 1 )
    beta_sim[ true_cov_beta ] <- true_coeff_beta*sample( c( 1 , -1 ), num_bal, replace = TRUE )
    
    y <- Balances%*%beta_sim + rnorm( subject_sim )
  }
  
  # Adjust inputs if x,y,z are provided
  Z_norm = z
  Z_norm[ Z_norm == 0 ] <- 0.5
  
  Z_norm <- Z_norm/rowSums( Z_norm )
  
  B_sim <- ifelse( is.null( z ), B_sim, nrow( z )  )
  covariates_sim <- ifelse( is.null( z ), covariates_sim, ncol( z ) )
  subject_sim <- ifelse( is.null( y ), subject_sim, dim( y )[1] )
  
  # Initiate starting values and allocate memory 
  samples <- floor( iterations/thin )
  
  # Intercept term alpha_j
  alpha. <- matrix( 0, nrow = B_sim, ncol = samples )
  
  # Inclusion indicators zeta_jp
  zeta. <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  
  # Regression Coefficients phi_jp
  phi. <- array( 0, dim = c( B_sim, covariates_sim, samples ) )
  
  # Psi for each i
  psi. <- array( 0, dim = c( subject_sim, B_sim, samples ) )
  
  # Inclusion indicators xi_m
  xi. <- matrix( 0, nrow = B_sim - 1, ncol = samples )
  
  # Adjust inital values for alpha, zeta, phi, psi, and xi if they are still NULL
  alpha.[ , 1] <- if( is.null( alpha ) ){ rnorm( n = 10 )}else{ alpha }  
  zeta.[ , , 1] <- ifelse( is.null( zeta ), 0, zeta )
  phi.[ , , 1]  <-  if( is.null( phi )){ rnorm( n = B_sim*covariates_sim )*zeta.[,,1]}else{ phi }
  psi.[ , , 1] <-  if( is.null( psi ) ){ as.matrix( Z_norm, nrow = subject_sim, ncol = B_sim) }else{ psi }
  xi.[ , 1] <- ifelse( is.null( xi ), 0, xi ) 
  
  # Run model
  output <- dm_lm_bvs( iterations, thin, alpha., y, z, x, phi., psi., sigma2_alpha, zeta., xi., sigma2_phi, a, b, a_0, b_0, h_alpha, h_beta, delta, rho, a_m, b_m)
  
  
  return( output )
}