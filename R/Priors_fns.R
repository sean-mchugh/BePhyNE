find_max_heights<-function(data_list, pred, bin_num){
  #datalist is a list of lists of your PA dtaa and occurrences over each species such that for one species
  #data_list[[1]]
  #       $species
  #       $y (1's and 0's for presences and absences)
  #       $ any number of subsequent entries for predictors, name them whatever but besure whatever they are named is what "pred" is defined as
  
  #ideally you do this after havign scaled your predictors, it will make your life a lot easier determining the bin size and may make it a lot more uniform between predictors
  #this only does a single bin size but you can put multiple predictors in, but if you need 2 predictors with different bin sizes I would say set it up two separate calls of this function and then make the outputs a list
  
  
  bin_size<-lapply(data_final, function(sp) diff(range(sp$X1[sp$y==1]))/bin_num)
  
  
  
  heights<-list()
  
  
  
  for(x in 1:length(pred)){
    
    bin_heights<-list()
    max_heights<-list()
    
    for(sp in 1:length(data_final)){
      
      
      bins<- seq(min(data_list[[sp]][[pred[[x]]]]), max(data_list[[sp]][[pred[[x]]]]), by = bin_size[[sp]])
      bin_heights[[sp]]<-list()
      
      for (i in 2:length(bins)){
        
        bin_heights[[sp]][[i]]<-sum(data_list[[sp]]$y[between(data_list[[sp]][[pred[[x]]]],bins[i-1], bins[i] )])/length(data_list[[sp]]$y[between(data_list[[sp]][[pred[[x]]]],bins[i-1], bins[i] )])
        
        
      }
      
      
      max_heights[[sp]]<-max(na.omit(unlist(bin_heights[[sp]])) )
      
    }
    
    
    
    heights[[x]]<-unlist(max_heights)
    
    #heights of actually 1 turn in inf when forward transformed
    
    heights[[x]][heights[[x]]==1]=.99
    
  }
  
  return(heights)
  
}




makePrior_ENE<- function(r, p, den.mu="unif", par.mu, den.sd="unif", par.sd, heights_by_sp=NULL, unif.corr=TRUE, Sigma=NULL, nu=NULL, plot=T){
  
  
  
  ## Save all parameters to use in subsequent functions.
  pars <- list()
  pars$r <- r ## Number of traits
  pars$p <- p ## Number of rate regimes.
  pars$den.mu <- den.mu
  pars$par.mu <-par.mu
  pars$den.sd <- den.sd
  pars$par.sd <- par.sd
  pars$unif.corr <- unif.corr
  pars$Sigma <- Sigma
  pars$nu <- nu
  
  
  par.mu.c <- par.mu[1,]
  par.mu.w <- par.mu[2,]
  
  pars$par.mu.c <- par.mu[1,]
  pars$par.mu.w <- par.mu[2,]
  
  
  if( nrow(par.mu)>2){
    
    par.mu.h <- par.mu[3,]
    pars$par.mu.h<- par.mu[3,]
    
  }
  
  if(plot==T){
    par(mfrow = c(2, 3))
  }
  
  
  ## Make a warning if 'unif.corr' is TRUE and 'Sigma' or 'nu' has values.
  if( unif.corr ){
    if( !is.null(Sigma) | !is.null(nu) ){
      warning("Found values for 'Sigma' and 'nu', but 'unif.corr=TRUE'. Using uniform correlations.")
    }
  }
  
  #trait pars
  #pars$c.min<-trait.ranges[1,][1]
  #pars$c.max<-trait.ranges[1,][2]
  #pars$w.min<-trait.ranges[2,][1]
  #pars$w.max<-trait.ranges[2,][2]
  #pars$h.min<-trait.ranges[3,][1]
  #pars$h.max<-trait.ranges[3,][2]
  
  ## Check if the distributions are known. This assumes all priors for the same parameter shares a distribution (with different parameters).
  ## Could extend this in the future.
  if(!den.mu %in% c("unif","norm") ) stop(" 'den.mu' need to be one of 'unif' or 'norm'.")
  if(!den.sd %in% c("unif","lnorm") ) stop(" 'den.sd' need to be one of 'unif' or 'lnorm'.")
  
  ## Check the number of parameters sent to the function. Those need to match.
  lpmu <- nrow(par.mu)
  if(!r == lpmu) stop("number of rows of 'par.mu' need to be equal to the number of traits in the model (r).")
  #if( is.vector(par.sd) && p > 1 ) stop("p need to be equal to the regimes fitted to the tree. So p == nrow(par.sd).")
  #if( is.matrix(par.sd) && !p == nrow(par.sd) ) stop("p need to be equal to the regimes fitted to the tree. So p == nrow(par.sd).")
  ## lpsd <- nrow(par.sd)
  ## if(!p == lpsd) stop("number of rows of 'par.sd' need to be equal to the number of R matrices fitted to the data (p).")
  if( unif.corr == FALSE ){
    if( is.matrix(Sigma) && p > 1 ) stop("Sigma need to be a list with elements equal to the number of regimes.")
    if( is.list(Sigma) && !p == length(Sigma) ) stop("length of 'Sigma' need to be equal to the number of regimes fitted to the tree. So length(Sigma) == p.")
    ## if( is.list(Sigma) && !p == length(Sigma) ) stop("length of 'Sigma' need to be equal to the length of 'nu'. So length(Sigma) == length(nu).")
    ## lsig <- length(Sigma)
    lnu <- length(nu)
    if(!p == lnu) stop("length of 'Sigma' and 'nu' need to be equal to the number of regimes fitted to the data. So length(Sigma) == length(nu) == p.")
  }
  
  ## Maybe add test for the correct parameter values? Would improve user experience.
  
  ## Note that the user will need to always specify the number of parameters equal to the number of the traits in the model.
  ## There was a problem in this part. The list of function with level 2 was having a strange behavior due to some environment thing.
  
  if(den.mu == "unif"){
    mn <- function(x) sum( sapply(1:length(x), function(i) stats::dunif(x[i], min=par.mu[i,1], max=par.mu[i,2], log=TRUE) ) )
  }
  if(den.mu == "norm"){
    
    if( is.null(heights_by_sp)==T){
      mn<- function(x) sum( c(stats::dnorm(x[1], mean=par.mu.c[[1]], sd=par.mu.c[[2]], log=TRUE),
                              stats::dlnorm(x[2],mean=par.mu.w[[1]], sd=par.mu.w[[2]], log=TRUE),
                              stats::dunif(x[3], min=par.mu.h[[1]], max=par.mu.h[[2]], log=TRUE)
      ) )
      
      if(plot==T){
        plot(density(stats::rnorm(10000, mean=par.mu[1,1], sd=par.mu[1,2]) ), main=paste("A_c") )
        plot(density(stats::rlnorm(10000, meanlog=par.mu[2,1], sdlog=par.mu[2,2]) ), main=paste("A_w")  )
        plot(density(stats::runif(10000, min=par.mu[3,1], max=par.mu[3,2]) ), main=paste("A_h")  )
        
      } }else{
        
        mn<- function(x) sum( c(stats::dnorm(x[1], mean=par.mu.c[[1]], sd=par.mu.c[[2]], log=TRUE),
                                stats::dlnorm(x[2],mean=par.mu.w[[1]], sd=par.mu.w[[2]], log=TRUE)
        ) )
        
        if(plot==T){
          plot(density(stats::rnorm(10000, mean=par.mu[1,1], sd=par.mu[1,2]) ), main=paste("A_c") )
          plot(density(stats::rlnorm(10000, meanlog=par.mu[2,1], sdlog=par.mu[2,2]) ), main=paste("A_w")  )
          
        }
      }
    
    
    
    
  }
  
  if( p == 1 ){
    if(den.sd == "unif"){
      sd <- function(x) sum(unlist(lapply(1:length(x), function(i)  stats::dunif(x[i], min=par.sd[[i,1]] , max=par.sd[[i,2]] , log=TRUE) )))
    }
    if(den.sd == "lnorm"){
      sd <- function(x) sum(unlist(lapply(1:length(x), function(i)  stats::dlnorm(x[i], meanlog=par.sd[[i,1]], sdlog=par.sd[i,2], log=TRUE) )))
      
      
      
      if(plot==T){
        plot(density(stats::rlnorm(10000, meanlog=par.sd[1,1], sdlog=par.sd[1,2]) ), main=paste("r_sd_c") )
        plot(density(stats::rlnorm(10000, meanlog=par.sd[2,1], sdlog=par.sd[2,2]) ), main=paste("r_sd_w") )
        if( is.null(heights_by_sp)==T){
          
          
          plot(density(stats::rlnorm(10000, meanlog=par.sd[3,1], sdlog=par.sd[3,2]) ), main=paste("r_sd_h")  )
          
        }
      }
      
      
    }
  } else{
    ## This thing below might solve the problem:
    if(den.sd == "unif"){
      ## The 'x' elements need to match with 'par.sd' rows.
      sd <- function(x) sum( sapply(1:p, function(i) sapply(x[[i]], function(y) stats::dunif(y, min=par.sd[i,1], max=par.sd[i,2], log=TRUE) ) ) )
    }
    if(den.sd == "lnorm"){
      ## The 'x' elements need to match with 'par.sd' rows.
      sd <- function(x) sum( sapply(1:p, function(i) sapply(x[[i]], function(y) stats::dlnorm(y, meanlog=par.sd[i,1], sdlog=par.sd[i,2], log=TRUE) ) ) )
      
    }
  }
  
  if(unif.corr == TRUE){
    if(p == 1){
      ## Only one matrix to work with.
      corr <- function(x) logDensityIWish(x, v=ncol(x)+1, diag(nrow=ncol(x)) )
    } else{
      ## When the prior on the correlation is set to uniform then all rate matrices will have the same prior.
      corr <- function(x) sum( sapply(x, function(y) logDensityIWish(y, v=ncol(y)+1, diag(nrow=ncol(y)) ) ) )
    }
  } else{
    if(p == 1){
      if(!is.matrix(Sigma)) stop(" 'Sigma' needs to be a matrix when a single regime is fitted to the data (p=1). ")
      if(!is.numeric(nu) || nu < r) stop(" 'nu' need to be a numeric value larger than the dimension of 'Sigma'. ")
      center <- (nu - ncol(Sigma) -1) * Sigma
      corr <- function(x) logDensityIWish(x, v=nu, center)
    } else{
      if(!all(sapply(Sigma, is.matrix))) stop(" 'Sigma' needs to be a list of matrices with length == p. ")
      if(!all(sapply(nu, is.numeric)) || all(nu < r) ) stop(" 'nu' needs to be a vector with numeric values larger than the dimension of R. ")
      corr_regime <- list()
      center <- list()
      for( i in 1:p){
        center[[i]] <- (nu[i] - ncol(Sigma[[i]]) -1) * Sigma[[i]]
        ## corr_regime[[i]] <- function(x) sapply(x, function(y) logDensityIWish(y, v=nu[i], center[[i]]) )
        corr_regime[[i]] <- function(x) logDensityIWish(x, v=nu[i], center[[i]])
      }
      corr <- function(x) sum( sapply(1:p, function(i) corr_regime[[i]](x[[i]]) ) )
    }
  }
  
  
  #just setting uninformative values based on range
  #expand to add other distributions
  #center.prior<- function(x) dunif(x, min=pars$c.min, max=pars$c.max)
  #width.prior<- function(x) dunif(x, min=pars$w.min, max=pars$w.max)
  #height.prior<-function(x) dunif(x, min=pars$h.min, max=pars$h.max)
  #sd<-function(x) dunif(x[1], 0, par.sd[[1]], log=T) + dunif(x[2], 0, par.sd[[2]], log=T) + dunif(x[3], 0, par.sd[[3]], log=T)
  
  ####height prior##
  #height is in the forward transformed space, allows us to use normal distribution for prior rather than beta dist.
  
  if(is.null(heights_by_sp)==T){
    
    height_priors=NULL
    
  }else{
    
    
    height_priors<-list()
    
    Y=((heights_by_sp-.05)/.95)
    FT_heights<- -1*log(Y/(1-Y))
    pars$FT_heights=FT_heights
    pars$heights=heights_by_sp
    
    #height_priors=list()
    #making a list of functions rather than one big function so we can call the specific function we want only when the move on that par happens,
    #do we really need to calculate the prior for EVERY sp EVERY iteration...nah thats gonna slow the whole mcmc down so instead lets just do it when there is a proposal on that height
    height_priors<-lapply(1:length(FT_heights), function(sp){ force(FT_heights[[sp]]);
      (function(x) dnorm(x, mean=FT_heights[[sp]], 0.15, log=T) )  })
    
    
  }
  
  
  
  
  
  
  
  
  
  res <- list(mean.prior=mn, corr.prior=corr, sd.prior=sd,  pars=pars, heights.prior=height_priors)
  class( res ) <- "prior_function"
  return( res )
}



make_all_priors <- function(
    N,
    tips,
    bd_range = NA,
    # ── μ (root) priors ──────────────────────────────────────────
    root_opt_mean         = rep(0,        N),
    root_opt_sd           = rep(0.2,      N),
    root_brdth_meanlog    = rep(log(.3),  N),
    root_brdth_sdlog      = rep(0.1,      N),
    # ── σ (log‑normal) hyper‑priors ──────────────────────────────
    sigsq_opt_meanlog     = rep(log(.1),  N),
    sigsq_opt_sdlog       = rep(0.5,      N),
    sigsq_brdth_meanlog   = rep(log(.1),  N),
    sigsq_brdth_sdlog     = rep(0.5,      N),
    # ── heights (default hard‑coded) ─────────────────────────────
    heights_by_sp         = sample(.95, size = tips, replace = TRUE),
    # ── constants forwarded to makePrior_ENE ─────────────────────
    r = 2, p = 1,
    plot = TRUE
){
  #users can specify a "breadth range that they expect species to cover for eahc predictor breadths ranges can be predictor spe cific or general
  if(sum(!is.na(bd_range))>0){
    if(class(bd_range)=="list"){
      for(i in length(bd_range)){
        max_bd<-max(bd_range[[i]])
        min_bd<-min(bd_range[[i]])
        bd<-((log(max_bd)-log(min_bd))/4)^2
        sigsq_brdth_meanlog[[i]]   = log(bd)
      }
    }else{
      max_bd<-max(bd_range)
      min_bd<-min(bd_range)
      bd<-((log(max_bd)-log(min_bd))/4)^2
      sigsq_brdth_meanlog   = rep(log(bd),  N)
      
    }
  }
  
  ## 0. Coerce heights into a list of length N ------------------
  height_list <- if (is.list(heights_by_sp)){
    stopifnot(length(heights_by_sp) == N)
    heights_by_sp
  } else {
    replicate(N, heights_by_sp, simplify = FALSE)
  }
  
  ## 1. Build priors for each character -------------------------
  prior_scale <- lapply(seq_len(N), function(i){
    
    ## μ matrix (2 × 2)
    par_mu <- matrix(c(root_opt_mean[i],      root_opt_sd[i],
                       root_brdth_meanlog[i], root_brdth_sdlog[i]),
                     nrow = 2, byrow = TRUE)
    
    ## σ matrix (2 × 2)
    par_sigsq <- matrix(c(sigsq_opt_meanlog[i],  sigsq_opt_sdlog[i],
                          sigsq_brdth_meanlog[i], sigsq_brdth_sdlog[i]),
                        nrow = 2, byrow = TRUE)
    
    makePrior_ENE(
      r   = r,  p = p,
      den.mu = "norm",
      heights_by_sp = height_list[[i]],
      par.mu = par_mu,
      den.sd = "lnorm",
      par.sd = par_sigsq,
      plot   = plot
    )
  })
  
  invisible(prior_scale)
}
