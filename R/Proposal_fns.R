#######proposal functions#########################################################################################################################################################################

make_tuning <- function(tree, pred,
                        center_slide = 0.18,
                        center_mult  = 0.12,
                        width_slide  = 0.15,
                        width_mult   = 0.2,
                        height_slide = 0.5,
                        height_mult  = 0.3,
                        w_mu  = c(0.6, 0.6),
                        w_sd  = c(0.15, 0.12),
                        v_cor = 100,
                        weights_height =2,
                        weights_center =3,
                        weights_width  =3,
                        weights_theta  =1,
                        weights_R_corr =1,
                        weights_R_sd   =1
) {

  n_species <- length(tree$tip.label)

  tuning <- list(
    niche_prop = lapply(1:pred, function(p) list(
      slide = tibble(
        center = sample(center_slide, n_species, replace = TRUE),
        width  = sample(width_slide,  n_species, replace = TRUE),
        height = sample(height_slide, n_species, replace = TRUE)
      ),
      multi = tibble(
        center = sample(center_mult, n_species, replace = TRUE),
        width  = sample(width_mult,  n_species, replace = TRUE),
        height = sample(height_mult, n_species, replace = TRUE)
      )
    )),
    w_mu = lapply(1:pred, function(p) list(
      slide = w_mu,
      multi = w_mu
    )),
    w_sd = lapply(1:pred, function(p) list(
      slide = w_sd,
      multi = w_sd
    )),
    v_cor = lapply(1:pred, function(p) v_cor)
  )

  move_weights= c(weights_height
                ,weights_center
                ,weights_width
                ,weights_theta
                ,weights_R_corr
                ,weights_R_sd  )

  move_prob=c("height" = move_weights[[1]]/sum(move_weights),
              "center" = move_weights[[2]]/sum(move_weights),
              "width"  = move_weights[[3]]/sum(move_weights),
              "theta"  = move_weights[[4]]/sum(move_weights),
              "R_corr" = move_weights[[5]]/sum(move_weights),
              "R_sd"   = move_weights[[6]]/sum(move_weights))



  return(list(tuning=tuning, move_probs=move_prob))
}



###niche prop###

propMVBM_by_clade <- function( tree, td, R, R_cor, R_sd, theta, n, curr_jac, V=NULL){
  if(is.null(V)){
    V <- ape::vcv.phylo(tree)
  }

  node =  sample((length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode  ) ,1)
  #extract clade to make shift
  clade<-extract.clade(tree, node = node )
  #td$dat[tree[[1]]$tip.label==clade$tip.label]
  ntraits <- ncol(td$dat)
  dat <- data.frame(species=tree$tip.label, td$dat)

  class(clade$tip.label[1])
  class(tree$tip.label[1])

  missing<-(1:nrow(dat))[tree$tip.label%in%clade$tip.label]

  old <- td$dat[missing,]
  dat[missing,] <- NA

  PPE <- Rphylopars::phylopars(trait_data = dat, tree = td$phy, phylocov_fixed = R, model_par_fixed = theta,
                               pheno_error = FALSE, skip_optim = TRUE, skip_EM = TRUE)
  ntraits <- ncol(td$dat)
  M <- matrix(0, nrow=ntraits*n, ncol=ntraits*n)
  for(i in 1:n){
    ind <- (1+ntraits*(i-1)):(ntraits*i)
    M[ind, ind] <- PPE$anc_cov[[missing[i]]]
  }
  cM <- chol(M)
  K <- t(cM) %*% kronecker(V[missing,missing], diag(1, nrow=ntraits, ncol=ntraits)) %*% cM
  A<- as.vector(t(PPE$anc_recon[missing,]))
  prop <- MASS::mvrnorm(1, A, K)
  td$dat[missing,] <- matrix(prop, nrow=n, ncol=ntraits, byrow = TRUE)
  hr <- mvtnorm::dmvnorm(prop, A, K, log=TRUE) - mvtnorm::dmvnorm(old, A, K, log=TRUE)
  return(list(td=td, hr=hr, missing=missing, n=n, K=K, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}

Slide_Proposal_byClade_byTrait_fn <- function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, d, niche_move){
  #d=.2
  # full_td<-td
  #td<-full_td
  #
  #tree=tree
  #td=td[[x]]
  #R=current_vals[[1]][[1]][[x]]$R
  #R_cor=current_vals[[1]][[1]][[x]]$R_cor
  #R_sd=current_vals[[1]][[1]][[x]]$R_sd
  #theta=current_vals[[1]][[1]][[x]]$A
  #n=n
  #curr_jac=current_vals[[1]][[1]][[x]]$curr.jac
  #H_fixed=H_fixed
  #d=d[[x]]
  node =  sample(1:(length(tree$tip.label)+tree$Nnode  ) ,1)
  #extract clade to make shift
  if(node <= length(tree$tip.label)){
    clade <- list(tip.label=tree$tip.label[node])
  } else {
    clade<-extract.clade(tree, node = node )
  }
  #td$dat[tree[[1]]$tip.label==clade$tip.label]
  ntraits <- ncol(td$dat)
  dat <- data.frame(species=tree$tip.label, td$dat)
  class(clade$tip.label[1])
  class(tree$tip.label[1])
  missing<-(1:nrow(dat))[tree$tip.label%in%clade$tip.label]
  old <- td$dat[missing,]
  if( niche_move=="center"){
    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,1]))
    hr <- 0
    td$dat[missing,1] <- matrix(prop, nrow=length(clade$tip.label), ncol=1, byrow = TRUE)
  } else if ( niche_move=="width"){
    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,2]))
    hr <- 0
    td$dat[missing,2] <- matrix(prop, nrow=length(clade$tip.label), ncol=1, byrow = TRUE)
  }else if ( niche_move=="height"){
    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,3]))
    hr <- 0
    td$dat[missing,3] <- matrix(prop, nrow=length(clade$tip.label), ncol=1, byrow = TRUE)
  }
  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}


Slide_Proposal_byClade_byTrait_fn <- function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, d, niche_move, clade){
  #d=.2
  # full_td<-td
  #td<-full_td
  #
  #tree=tree
  #td=td[[1]]
  #R=current_vals[[1]][[1]][[x]]$R
  #R_cor=current_vals[[1]][[1]][[x]]$R_cor
  #R_sd=current_vals[[1]][[1]][[x]]$R_sd
  #theta=current_vals[[1]][[1]][[x]]$A
  #n=n
  #curr_jac=current_vals[[1]][[1]][[x]]$curr.jac
  #H_fixed=H_fixed
  #d=d[[x]]
  #node =  sample(1:(length(tree$tip.label)+tree$Nnode  ) ,1)
  ##extract clade to make shift
  #if(node <= length(tree$tip.label)){
  #  clade <- list(tip.label=tree$tip.label[node])
  #} else {
  #  clade<-extract.clade(tree, node = node )
  #}
  ##td$dat[tree[[1]]$tip.label==clade$tip.label]

  ntraits <- ncol(td$dat)

  dat <- data.frame(species=tree$tip.label, td$dat)
  #class(clade$tip.label[1])
  #class(tree$tip.label[1])
  missing<-(1:nrow(dat))[tree$tip.label%in%clade$tip.label]
  old <- td$dat[missing,]



  # d[[1]]$slide[missing,1]

  if( niche_move=="center"){
    prop <-  d * (stats::runif(1) - 0.5) + as.vector(t(old[,1]))
    hr <- 0
    td$dat[missing,1] <-  prop

  } else if ( niche_move=="width"){
    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,2]))
    hr <- 0
    td$dat[missing,2] <-prop
  }else if ( niche_move=="height"){
    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,3]))
    hr <- 0
    td$dat[missing,3] <- prop
  }
  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}


Slide_Proposal_fn <- function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, d, niche_move){
  #d=.2
  ntraits <- ncol(td$dat)
  dat <- data.frame(species=tree$tip.label, td$dat)
  missing <- sample(1:nrow(dat), n, replace=FALSE)
  old <- td$dat[missing,]

  if( niche_move=="center"){

    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,1]))
    hr <- 0

    td$dat[missing,1] <- matrix(prop, nrow=n, ncol=1, byrow = TRUE)


  } else if ( niche_move=="width"){

    prop <- d * (stats::runif(1) - 0.5) + as.vector(t(old[,2]))
    hr <- 0

    td$dat[missing,2] <- matrix(prop, nrow=n, ncol=1, byrow = TRUE)


  }
  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}


Multiplier_Proposal_byClade_byTrait_fn <- function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, d, niche_move){
  # d=.2

  node =  sample((length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode  ) ,1)
  #extract clade to make shift
  clade<-extract.clade(tree, node = node )
  #td$dat[tree[[1]]$tip.label==clade$tip.label]

  dat <- data.frame(species=tree$tip.label, td$dat)
  missing<-(1:nrow(td$dat))[tree$tip.label%in%clade$tip.label]

  old <- td$dat[missing,]
  m <- exp(d * (stats::runif(1) - 0.5))
  #prop <- as.vector(t(old)) * m
  #hr <- log(m)
  if( niche_move=="center"){

    prop <- as.vector(t(old[,1])) * m
    hr <- log(m)

    td$dat[missing,1] <- matrix(prop, nrow=length(clade$tip.label), ncol=1, byrow = TRUE)

  } else if( niche_move=="width") {

    prop <- as.vector(t(old[,2])) * m
    hr <- log(m)

    td$dat[missing,2] <- matrix(prop, nrow=length(clade$tip.label), ncol=1, byrow = TRUE)

  }else if( niche_move=="height") {

    prop <- as.vector(t(old[,3])) * m
    hr <- log(m)

    td$dat[missing,3] <- matrix(prop, nrow=length(clade$tip.label), ncol=1, byrow = TRUE)

  }

  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}



Multiplier_Proposal_byClade_byTrait_fn <- function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, d, niche_move, clade){
  # d=.2

  #node =  sample((length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode  ) ,1)
  ##extract clade to make shift
  #clade<-extract.clade(tree, node = node )
  ##td$dat[tree[[1]]$tip.label==clade$tip.label]

  dat <- data.frame(species=tree$tip.label, td$dat)
  missing<-(1:nrow(td$dat))[tree$tip.label%in%clade$tip.label]

  old <- td$dat[missing,]

  m <- exp(d * (stats::runif(1) - 0.5))

  # m <- exp(d[[x]]$multi[missing,1] * (stats::runif(1) - 0.5))



  #prop <- as.vector(t(old)) * m
  #hr <- log(m)
  if( niche_move=="center"){

    prop <- as.vector(t(old[,1])) * m
    hr <- log(m)

    td$dat[missing,1] <- prop

  } else if( niche_move=="width") {

    prop <- as.vector(t(old[,2])) * m
    hr <- log(m)

    td$dat[missing,2] <- prop

  }else if( niche_move=="height") {

    prop <- as.vector(t(old[,3])) * m
    hr <- log(m)

    td$dat[missing,3] <- prop

  }

  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}



Multiplier_Proposal_fn <- function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, d, niche_move){
  # d=.2
  ntraits <- ncol(td$dat)
  dat <- data.frame(species=tree$tip.label, td$dat)
  missing <- sample(1:nrow(dat), n, replace=FALSE)
  old <- td$dat[missing,]
  m <- exp(d * (stats::runif(1) - 0.5))
  #prop <- as.vector(t(old)) * m
  #hr <- log(m)
  if( niche_move=="center"){

    prop <- as.vector(t(old[,1])) * m
    hr <- log(m)

    td$dat[missing,1] <- matrix(prop, nrow=n, ncol=1, byrow = TRUE)

  } else if( niche_move=="width") {

    prop <- as.vector(t(old[,2])) * m
    hr <- log(m)

    td$dat[missing,2] <- matrix(prop, nrow=n, ncol=1, byrow = TRUE)

  }

  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}


###Evo Model prop###

Posdef <- function (n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

riwish <- function(v, S){
  S <- solve(S)
  if (!is.matrix(S)) S <- matrix(S)
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(stats::rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- stats::rnorm(p *(p - 1)/2)
  }
  out <- crossprod(Z %*% CC)
  return(solve(out))
}

makePropIWish_C <- function(vcv, k, v) {
  .Call('_ratematrix_makePropIWish_C', PACKAGE = 'ratematrix', vcv, k, v)
}

#hotfix for .call being pain on HPC, use for HPCs
makePropIWish_C <- function(curr.vcv, k, v){
  ## Make a proposal for a new vcv matrix given the current state.
  ## k = dimension of the matrix.
  ## v = degrees of freedom.
  center <- (v-k-1) * curr.vcv
  prop <- riwish(v, center)
  return(prop)
}


logDensityIWish <- function(W, v, S){
  ## This function is derived from the 'diwish' from MCMCpack.
  ## This returns the log density.
  k <- nrow(S)
  lgammapart <- 0
  for (i in 1:k) {
    lgammapart <- lgammapart + lgamma((v + 1 - i)/2)
  }
  ldenom <- lgammapart + ( ( (v*k)/2 ) * log( 2 ) ) + ( ( (k*(k-1))/4 ) * log( pi ) )
  detS <- det(S)
  detW <- det(W)
  hold <- S %*% solve(W)
  tracehold <- sum(hold[row(hold) == col(hold)])
  lnum <- ( ( v/2 ) * log( detS ) ) + ( ( -(v + k + 1)/2 ) * log( detW ) ) + ( -1/2 * tracehold )
  return(lnum - ldenom)
}

cov2cor_C <- function(V) {
  .Call('_ratematrix_cov2cor_C', PACKAGE = 'ratematrix', V)
}

##hotfix for .call being pain on HPC, use for HPCs

cov2cor_C<-function(V){
  return(cov2cor(V))
}

slideWindow <- function(x, w){
  ## Sliding window proposal for unbounded trait.
  ## x = the current value.
  ## w = the width parameter of the proposal.
  y <- stats::runif(1, min = x - (w/2), max = x + (w/2) )
  return(y)
}

multiplierProposal <- function(x, a){
  ## This proposal scheme will perform steps in the log space of the parameter.
  ## This is much more efficient for the cases in which small changes in large
  ##    values effect less than small changes in small values.
  ## Note this proposal will never change the sign of the parameter.
  ## x = the current value.
  ## a = the scale of the multiplier.
  ## lambda <- 2 * log(a)
  ## m <- exp( lambda * runif(1, min=-0.5, max=0.5) )

  ## The previous multiplier proposal was assuming the parameter for the model was in log space.
  ## This is the correct implementation for the multiplier proposal.
  m <- exp( a * (runif(1) - 0.5) )
  ## Note that 'm' here is the proposal ratio. So need to spit this out.
  return( setNames(c(m * x, m), c("prop","prop.ratio") ) )
}

slideWindow_theta <- function(w_mu, tree, td, R, R_cor, R_sd, theta, n, curr_jac){
  ## Sliding window proposal for unbounded trait.
  ## x = the current value.
  ## w = the width parameter of the proposal.
  prop.theta <- sapply(1:length(theta), function(x) slideWindow(theta[x], w_mu[x]) )

  return(list(td=td, hr=0, missing=missing, n=n, K=NULL, A=prop.theta,  R=R, R_cor=R_cor, R_sd=R_sd, prop_jac=curr_jac, jj=0))
}

SD_prop<-function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, w_sd){
  num_error=0
  repeat {
    num_error <- num_error+1


    multi.sd <- sapply(1:length(R_sd), function(x) multiplierProposal(x = R_sd[[x]], a = w_sd[[x]]) )

    prop.sd <- multi.sd[1,]

    prop.R <- rebuild.cov( r=R_cor, v=prop.sd^2)

    if (num_error==10){
      return(NULL)
    }

    if (is.positive.definite(prop.R )==T){
      break
    }
  }
  hr=sum(log(multi.sd[2,]))

  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=prop.R, R_cor=R_cor, R_sd=prop.sd, prop_jac=curr_jac,  jj=0))
}

SD_prop_slide<-function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, w_sd){
  num_error=0
  repeat {
    num_error <- num_error+1


    multi.sd <- sapply(1:length(R_sd), function(x) slideWindow(R_sd[[x]],  w_sd[[x]]) )

    prop.sd <- multi.sd

    prop.R <- rebuild.cov( r=R_cor, v=prop.sd^2)

    if (num_error==10){
      return(NULL)
    }

    if (is.positive.definite(prop.R )==T){
      break
    }
  }
  hr=0
  return(list(td=td, hr=hr, missing=missing, n=n, K=NULL, A=theta, R=prop.R, R_cor=R_cor, R_sd=prop.sd, prop_jac=curr_jac,  jj=0))
}



Cor_matrix_prop<-function(tree, td, R, R_cor, R_sd, theta, n, curr_jac, v=100){
  #may have been the big error 5/9 R_sd was used for jacobnian instead of variances
  #had make.symmetric around the rebuild, dont have that elsewhere, trying to see if removing it improves the MCMC add back if things come out worse maybe
  num_error=0

  repeat {

    num_error <- num_error+1
    prop_IW<-makePropIWish_C(R, length(2), v)
    prop_cor=cov2cor_C(prop_IW)
    prop.R <- rebuild.cov( r=prop_cor, v=R_sd^2)

    if (num_error==10){
      return(NULL)
    }

    if (is.positive.definite(prop.R )==T){
      break
    }
  }


  hr <- hastingsDensity(curr.vcv=R, prop.vcv=prop.R, p=length(R_sd), v=v)
  prop.r.jacobian <- sum( sapply(1:length(R_sd), function(x) log( ((R_sd[x])^2)) ) ) * log( (length(R_sd)-1)/2 )
  jj <- prop.r.jacobian - curr_jac
  return(list(td=td, hr=hr, missing=NULL, n=n, K=NULL, A=theta, R=prop.R, R_cor=prop_cor, R_sd=R_sd, prop_jac=prop.r.jacobian, jj=jj))
}


make.symmetric <- function(s){
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  s
} #helps fix some of the numerical evaluation results


####full prop fn#######

prop_fn<-function(prop, td, current_vals, n, d) {
  if (prop==1){
    ##center moves
    niche_move=sample(c(1,2,3), size=1, prob = c(0,1/2,1/2))



    if(niche_move==1){
      #proposal = lapply(1:length(td), function(x) propMVBM(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, V=NULL))
      move="Center_BM"

    } else if(niche_move==2){
      #proposal= lapply(1:length(td), function(x) Slide_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center, niche_move = "center"))
      proposal= lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center_slide, niche_move = "center"))
      if(center_fixed!=FALSE){
        for (pred in center_fixed){
          proposal[[pred]]$td<-td[[pred]]
        }
      }

      #td[[x]]$dat
      #proposal[[x]]$td$dat
      #proposal[[x]]$missing

      move="Center_Slide"

    }else if(niche_move==3){
      #proposal = lapply(1:length(td), function(x) Multiplier_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center, niche_move = "center"))
      proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center_mult, niche_move = "center"))
      if(center_fixed!=FALSE){
        for (pred in center_fixed){
          proposal[[pred]]$td<-td[[pred]]
        }
      }

      move="Center_Multiplier"
      #lol<-unlist(proposal, recursive=F)
    }
  }else if (prop==2){
    ##width moves
    niche_move=sample(c(1,2,3), size=1, prob = c(0,1/2,1/2))


    if(niche_move==1){
      #proposal = lapply(1:length(td), function(x) propMVBM(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, V=NULL))
      move="Width_BM"

    } else if(niche_move==2){
      #proposal= lapply(1:length(td), function(x) Slide_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, H_fixed = H_fixed, d[[x]]))
      proposal= lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$width_slide, niche_move="width"))

      # td[[x]]$dat== proposal[[x]]$td$dat
      #proposal[[x]]$td$dat
      move="Width_Slide"

    }else if(niche_move==3){
      #proposal = lapply(1:length(td), function(x) Multiplier_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, H_fixed = H_fixed, d[[x]]))
      proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$width_mult, niche_move="width"))

      move="Width_Multiplier"
      #lol<-unlist(proposal, recursive=F)
    }


  } else if(prop==3){
    ##theta move
    proposal <- lapply(1:length(td), function(x) slideWindow_theta(w_mu[[x]], tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac))
    move="theta"

  } else if(prop==4){
    ##correlation matrix component of R matrix move
    proposal <- lapply(1:length(td), function(x) Cor_matrix_prop(tree,td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, v=v_cor[[x]]))
    move="R_corr"

    if (is.null(proposal)==T){
      print("Your matrix is fucked yo, this is the signularity error")
      if(trim==TRUE){

        return(list(acceptances=acceptances, accepted.moves=list(A.dat.bm.moves=A.dat.bm.moves, A.dat.slide.moves=A.dat.slide.moves, A.dat.mult.moves=A.dat.mult.moves, A.theta.moves=A.theta.moves, A.R.corr.moves= A.R.corr.moves, A.R.sd.moves=A.R.sd.moves), rejected.moves=list(Rej.dat.bm.moves=Rej.dat.bm.moves, Rej.dat.slide.moves=Rej.dat.slide.moves, Rej.dat.mult.moves=Rej.dat.mult.moves, Rej.theta.moves=Rej.theta.moves, Rej.R.corr.moves= Rej.R.corr.moves, Rej.R.sd.moves=Rej.R.sd.moves), trimmed_chain=trimmed_chain, NAs=NA_moves, NA.moves=list(NA.dat.bm.moves=NA.dat.bm.moves, NA.dat.slide.moves=NA.dat.slide.moves, NA.dat.mult.moves=NA.dat.mult.moves, NA.theta.moves=NA.theta.moves, NA.R.corr.moves= NA.R.corr.moves, NA.R.sd.moves=NA.R.sd.moves), "full"))

      }else{

        return(list(acceptances=acceptances, accepted.moves=list(A.dat.bm.moves=A.dat.bm.moves, A.dat.slide.moves=A.dat.slide.moves, A.dat.mult.moves=A.dat.mult.moves, A.theta.moves=A.theta.moves, A.R.corr.moves= A.R.corr.moves, A.R.sd.moves=A.R.sd.moves), rejected.moves=list(Rej.dat.bm.moves=Rej.dat.bm.moves, Rej.dat.slide.moves=Rej.dat.slide.moves, Rej.dat.mult.moves=Rej.dat.mult.moves, Rej.theta.moves=Rej.theta.moves, Rej.R.corr.moves= Rej.R.corr.moves, Rej.R.sd.moves=Rej.R.sd.moves), chain=chain, NAs=NA_moves, NA.moves=list(A.dat.bm.moves=NA.dat.bm.moves, NA.dat.slide.moves=NA.dat.slide.moves, NA.dat.mult.moves=NA.dat.mult.moves, NA.theta.moves=NA.theta.moves, NA.R.corr.moves= NA.R.corr.moves, NA.R.sd.moves=NA.R.sd.moves), "full"))

      }
    }

  } else if(prop==5){
    ##vector of standard deviations component of R matrix move
    proposal <- lapply(1:length(td), function(x) SD_prop(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, w_sd=w_sd[[x]]))
    move="R_sd"
    if (is.null(proposal)==T){
      print("Your matrix is fucked yo, this is the signularity error")
      if(trim==TRUE){
        return(list(acceptances=acceptances, accepted.moves=list(A.dat.bm.moves=A.dat.bm.moves, A.dat.slide.moves=A.dat.slide.moves, A.dat.mult.moves=A.dat.mult.moves, A.theta.moves=A.theta.moves, A.R.corr.moves= A.R.corr.moves, A.R.sd.moves=A.R.sd.moves), rejected.moves=list(Rej.dat.bm.moves=Rej.dat.bm.moves, Rej.dat.slide.moves=Rej.dat.slide.moves, Rej.dat.mult.moves=Rej.dat.mult.moves, Rej.theta.moves=Rej.theta.moves, Rej.R.corr.moves= Rej.R.corr.moves, Rej.R.sd.moves=Rej.R.sd.moves), trimmed_chain=trimmed_chain, NAs=NA_moves, NA.moves=list(NA.dat.bm.moves=NA.dat.bm.moves, NA.dat.slide.moves=NA.dat.slide.moves, NA.dat.mult.moves=NA.dat.mult.moves, NA.theta.moves=NA.theta.moves, NA.R.corr.moves= NA.R.corr.moves, NA.R.sd.moves=NA.R.sd.moves), "full"))
      }else{
        return(list(acceptances=acceptances, accepted.moves=list(A.dat.bm.moves=A.dat.bm.moves, A.dat.slide.moves=A.dat.slide.moves, A.dat.mult.moves=A.dat.mult.moves, A.theta.moves=A.theta.moves, A.R.corr.moves= A.R.corr.moves, A.R.sd.moves=A.R.sd.moves), rejected.moves=list(Rej.dat.bm.moves=Rej.dat.bm.moves, Rej.dat.slide.moves=Rej.dat.slide.moves, Rej.dat.mult.moves=Rej.dat.mult.moves, Rej.theta.moves=Rej.theta.moves, Rej.R.corr.moves= Rej.R.corr.moves, Rej.R.sd.moves=Rej.R.sd.moves), chain=chain, NAs=NA_moves, NA.moves=list(A.dat.bm.moves=NA.dat.bm.moves, NA.dat.slide.moves=NA.dat.slide.moves, NA.dat.mult.moves=NA.dat.mult.moves, NA.theta.moves=NA.theta.moves, NA.R.corr.moves= NA.R.corr.moves, NA.R.sd.moves=NA.R.sd.moves), "full"))

      }
    }
  }


  return(proposal)

}

