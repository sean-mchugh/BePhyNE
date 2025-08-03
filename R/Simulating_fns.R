#simulation_fn#####################################################################################################################################################################

Posdef <- function (n, ev = runif(n, 0, 10))
{
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


priorSim_pars<-function(Prior, phylo, dist, heights){
  #R_true[[i]] <- matrix(c(10, 0,  0,
  #                        0,  .1,  0,
  #                        0,  0, .1), byrow=TRUE,ncol=3, nrow=3)
  #R_decom<-decompose.cov(R_true[[i]])
  #R_corr_true[[i]]<-R_decom$r
  #R_sd_true[[i]]<-sqrt(R_decom$v)
  if(dist=="unif"){
    
    R_sd_true<-lapply(1:length(Prior), function(i) c(runif(1, Prior[[i]]$pars$par.sd[1,1], Prior[[i]]$pars$par.sd[[1,2]]), runif(1, Prior[[i]]$pars$par.sd[2,1], Prior[[i]]$pars$par.sd[2,2])))
    
    
    
    R_cor_true_R<-lapply(1:length(Prior), function(i) riwish(3,diag(2)))
    
    R_cor_true=lapply(1:length(Prior), function(i) cov2cor_C(R_cor_true_R[[i]]))
    
    R_true<-lapply(1:length(Prior), function(i) rebuild.cov(R_cor_true[[i]],(R_sd_true[[i]]^2)))
    
    
    # A_true[[i]]= c(20,3,.8)
    #tA[[i]]= forwardTransform1(A_true[[i]])
    #[1] 14.000000  1.609438 -1.321756
    
    A_true_full=lapply(1:length(Prior), function(i) forwardTransform1(c(runif(1,
                                                                              min = Prior[[i]]$pars$par.mu[1,1],
                                                                              max = Prior[[i]]$pars$par.mu[1,2]),
                                                                        runif(1,
                                                                              min = Prior[[i]]$pars$par.mu[2,1],
                                                                              max = Prior[[i]]$pars$par.mu[2,2]))
    ))
  } else if (dist=="norm"){
    
    R_sd_true<-lapply(1:length(Prior), function(i) c(rlnorm(n=1, meanlog=Prior[[i]]$pars$par.sd[1,1], sdlog=Prior[[i]]$pars$par.sd[[1,2]]), rlnorm(1, Prior[[i]]$pars$par.sd[2,1], Prior[[i]]$pars$par.sd[2,2])))
    
    R_cor_true_R<-lapply(1:length(Prior), function(i) riwish(3,diag(2)))
    
    R_cor_true=lapply(1:length(Prior), function(i) cov2cor_C(R_cor_true_R[[i]]))
    
    R_true<-lapply(1:length(Prior), function(i) rebuild.cov(R_cor_true[[i]],(R_sd_true[[i]]^2)))
    
    
    # A_true[[i]]= c(20,3,.8)
    #tA[[i]]= forwardTransform1(A_true[[i]])
    #[1] 14.000000  1.609438 -1.321756
    
    A_true_full=lapply(1:length(Prior), function(i) forwardTransform1(c(rnorm(1,
                                                                              Prior[[i]]$pars$par.mu[1,1], # mean
                                                                              Prior[[i]]$pars$par.mu[1,2]),# sd
                                                                        rlnorm(1,
                                                                               Prior[[i]]$pars$par.mu[2,1],  # mean
                                                                               Prior[[i]]$pars$par.mu[2,2])) # sd
    ))
    
  }
  
  
  
  A_true_bt_full<-lapply( 1:length(Prior), function(i) backTransform1(A_true_full[[i]]))
  
  
  A_true<-lapply(1:length(Prior), function(i) A_true_full[[i]][1:2])
  A_true_bt<-lapply(1:length(Prior), function(i) A_true_bt_full[[i]][1:2])
  
  true_dat_no_H <- lapply(1:length(Prior), function(i) mvMORPH::mvSIM(phylo, nsim = 1, param=list(ntraits=2, sigma=R_true[[i]], theta=A_true[[i]])))
  sim_dat<- lapply(1:length(Prior), function(i) cbind(true_dat_no_H[[i]], as.numeric(seq(-2, -2, length.out=length(phylo$tip.label)))))
  
  
  
  for ( pred in 1:length(Prior)){
    
    
    Y=((heights[[pred]]-.05)/.95)
    FT_heights<- -1*log(Y/(1-Y))
    
    sim_dat[[pred]][,3]<-FT_heights
  }
  
  
  sim_dat_bt=lapply(1:length(sim_dat), function(x) t(apply(sim_dat[[x]], 1, backTransform1)))
  
  sim_td<-lapply(1:length(sim_dat), function(x) make.treedata(phylo, sim_dat[[x]])$dat)
  
  
  sim_td_bt<-lapply(1:length(sim_dat_bt), function(x) make.treedata(phylo, sim_dat_bt[[x]])$dat)
  
  
  
  
  
  return(list( sim_dat=list(sim_dat_ft=sim_dat,
                            sim_dat_bt=sim_dat_bt,
                            sim_td=sim_td,
                            sim_td_bt=sim_td_bt),
               
               A=list(A_ft=A_true,
                      A_bt=A_true_bt),
               
               R=list(R=R_true,
                      R_sd=R_sd_true,
                      R_cor=R_cor_true )))
  
  
  
}

priorSim_pars<-function(Prior, phylo, dist, hard_coded_heights=NULL){
  #R_true[[i]] <- matrix(c(10, 0,  0,
  #                        0,  .1,  0,
  #                        0,  0, .1), byrow=TRUE,ncol=3, nrow=3)
  #R_decom<-decompose.cov(R_true[[i]])
  #R_corr_true[[i]]<-R_decom$r
  #R_sd_true[[i]]<-sqrt(R_decom$v)
  if(dist=="unif"){
    
    R_sd_true<-lapply(1:length(Prior), function(i) c(runif(1, Prior[[i]]$pars$par.sd[1,1], Prior[[i]]$pars$par.sd[[1,2]]), runif(1, Prior[[i]]$pars$par.sd[2,1], Prior[[i]]$pars$par.sd[2,2])))
    
    
    
    R_cor_true_R<-lapply(1:length(Prior), function(i) riwish(3,diag(2)))
    
    R_cor_true=lapply(1:length(Prior), function(i) cov2cor_C(R_cor_true_R[[i]]))
    
    R_true<-lapply(1:length(Prior), function(i) rebuild.cov(R_cor_true[[i]],(R_sd_true[[i]]^2)))
    
    
    # A_true[[i]]= c(20,3,.8)
    #tA[[i]]= forwardTransform1(A_true[[i]])
    #[1] 14.000000  1.609438 -1.321756
    
    A_true_full=lapply(1:length(Prior), function(i) forwardTransform1(c(runif(1,
                                                                              min = Prior[[i]]$pars$par.mu[1,1],
                                                                              max = Prior[[i]]$pars$par.mu[1,2]),
                                                                        runif(1,
                                                                              min = Prior[[i]]$pars$par.mu[2,1],
                                                                              max = Prior[[i]]$pars$par.mu[2,2]))
    ))
  } else if (dist=="norm"){
    
    R_sd_true<-lapply(1:length(Prior), function(i) unlist( lapply(1:nrow(Prior[[i]]$pars$par.sd), function(trait) c(rlnorm(n=1, meanlog=Prior[[i]]$pars$par.sd[trait,1], sdlog=Prior[[i]]$pars$par.sd[[trait,2]])))))
    
    R_cor_true_R<-lapply(1:length(Prior), function(i) riwish(nrow(Prior[[i]]$pars$par.sd)+1,diag(nrow(Prior[[i]]$pars$par.sd))))
    
    R_cor_true=lapply(1:length(Prior), function(i) cov2cor_C(R_cor_true_R[[i]]))
    
    R_true<-lapply(1:length(Prior), function(i) rebuild.cov(R_cor_true[[i]],(R_sd_true[[i]]^2)))
    
    
    # A_true[[i]]= c(20,3,.8)
    #tA[[i]]= forwardTransform1(A_true[[i]])
    #[1] 14.000000  1.609438 -1.321756
    
    if( is.null(hard_coded_heights)==T){
      
      
      A_true_full=lapply(1:length(Prior), function(i) forwardTransform1(c(rnorm(1,
                                                                                Prior[[i]]$pars$par.mu[1,1], # mean
                                                                                Prior[[i]]$pars$par.mu[1,2]),# sd
                                                                          rlnorm(1,
                                                                                 Prior[[i]]$pars$par.mu[2,1],  # mean
                                                                                 Prior[[i]]$pars$par.mu[2,2]),
                                                                          runif(1,
                                                                                Prior[[i]]$pars$par.mu[3,1],  # mean
                                                                                Prior[[i]]$pars$par.mu[3,2])) # sd
                                                                        # sd
      ))
      
      
    } else{
      
      A_true_full=lapply(1:length(Prior), function(i) forwardTransform1(c(rnorm(1,
                                                                                Prior[[i]]$pars$par.mu[1,1], # mean
                                                                                Prior[[i]]$pars$par.mu[1,2]),# sd
                                                                          rlnorm(1,
                                                                                 Prior[[i]]$pars$par.mu[2,1],  # mean
                                                                                 Prior[[i]]$pars$par.mu[2,2]))
      ))
      
      
      
    }
    
  }
  
  
  
  #check to see if you included ehight or not in mvBM, should expand to allow any comination of center,width, and height
  
  A_true_bt_full<-lapply( 1:length(Prior), function(i) backTransform1(A_true_full[[i]]))
  
  
  if(nrow(R_true[[1]])==2){
    
    A_true<-lapply(1:length(Prior), function(i) A_true_full[[i]][1:2])
    A_true_bt<-lapply(1:length(Prior), function(i) A_true_bt_full[[i]][1:2])
    
    true_dat_no_H <- lapply(1:length(Prior), function(i) mvMORPH::mvSIM(phylo, nsim = 1, param=list(ntraits=2, sigma=R_true[[i]][1:2,1:2], theta = A_true[[i]][1:2])))
    sim_dat<- lapply(1:length(Prior), function(i) cbind(true_dat_no_H[[i]], as.numeric(seq(-2, -2, length.out=length(phylo$tip.label)))))
    
    if(is.null(hard_coded_heights)==F){
      
      for ( pred in 1:length(Prior)){
        
        
        Y=((hard_coded_heights[[pred]]-.05)/.95)
        FT_hard_coded_heights<- -1*log(Y/(1-Y))
        
        sim_dat[[pred]][,3]<-FT_hard_coded_heights
      }
      
    }
    
    
  } else if(nrow(R_true[[1]]==3)){
    A_true<-lapply(1:length(Prior), function(i) A_true_full[[i]][1:3])
    A_true_bt<-lapply(1:length(Prior), function(i) A_true_bt_full[[i]][1:3])
    
    true_dat<- lapply(1:length(Prior), function(i) mvMORPH::mvSIM(phylo, nsim = 1, param=list(ntraits=nrow(R_true[[i]]), sigma=R_true[[i]], theta=A_true[[i]])))
    sim_dat<- lapply(1:length(Prior), function(i) cbind(true_dat[[i]]))
    
  }
  
  
  
  sim_dat_bt=lapply(1:length(sim_dat), function(x) t(apply(sim_dat[[x]], 1, backTransform1)))
  
  sim_td<-lapply(1:length(sim_dat), function(x) make.treedata(phylo, sim_dat[[x]])$dat)
  
  
  sim_td_bt<-lapply(1:length(sim_dat_bt), function(x) make.treedata(phylo, sim_dat_bt[[x]])$dat)
  
  
  
  
  
  return(list( sim_dat=list(sim_dat_ft=sim_dat,
                            sim_dat_bt=sim_dat_bt,
                            sim_td=sim_td,
                            sim_td_bt=sim_td_bt),
               
               A=list(A_ft=A_true,
                      A_bt=A_true_bt),
               
               R=list(R=R_true,
                      R_sd=R_sd_true,
                      R_cor=R_cor_true )))
  
  
  
}



findBadStart<- function(res, pa_data, plot=F){
  
  betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
  
  yy = lapply( 1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1+ betas[[1]][i,3]*(pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2+ betas[[2]][i,3]*(pa_data[[i]]$X2^2) )))
  
  presProb <-lapply(1:nrow(res[[1]]), function(i) 1/(1+exp(-1*yy[[i]])) )
  
  likelihood= lapply(1:nrow(res[[1]]), function(z) dbinom(pa_data[[z]]$y,1,prob=presProb[[z]],log=T))
  sumlikelihood=lapply(1:length(likelihood), function(x)  sum(likelihood[[x]]))
  
  ##identify which species have -inf starting params
  inf<-(1:length(sumlikelihood))[sumlikelihood==-Inf]
  
  return(inf)
  
}


findStart<-function(res, pa_data, tree, plot=F){
  
  #set.seed(1)
  
  #often times when simulating data/ you get underflow and a flat likelihood floor, to get start values that dont hit the floor we are brute forcing
  #the center proposal iteratively for each species until you have working starting parameters
  #this may need to be used also if trying to use ML start params as optim would have the same issue as our estimator
  res_new<-res
  
  betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
  
  yy = lapply( 1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1+ betas[[1]][i,3]*(pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2+ betas[[2]][i,3]*(pa_data[[i]]$X2^2) )))
  
  presProb <-lapply(1:nrow(res[[1]]), function(i) 1/(1+exp(-1*yy[[i]])) )
  
  likelihood= lapply(1:nrow(res[[1]]), function(z) return(dbinom(pa_data[[z]]$y,1,prob=presProb[[z]],log=T)))
  sumlikelihood=lapply(1:length(likelihood), function(x)  sum(likelihood[[x]]))
  
  ##identify which species have -inf starting params
  inf<-(1:length(sumlikelihood))[sumlikelihood==-Inf]
  
  if(length(inf)==0){
    #paste("no starting pars in the flat lands")
    res_ft<-lapply(1:length(res), function(x) t(apply(res[[x]], 1, forwardTransform1)) )
    res_td<-lapply(1:length(res), function(x) make.treedata(tree, res_ft[[x]])$dat)
    
    res_td_bt<-lapply(1:length(res), function(x) make.treedata(tree, res[[x]])$dat )
    
    return(list(start_pars_bt=res, start_pars_ft=res_ft, start_td=res_td, start_td_bt=res_td_bt, likelihoods=sumlikelihood))
  }
  
  sp <- lapply(1:length(res), function(x) res[[x]][inf,] )
  sp_old<-sp
  
  pa_data_inf<-pa_data[inf]
  #View(pa_data_inf)
  
  
  
  sum=list()
  
  ##iterate over each species
  for (y in 1:length(inf)){
    ##count to keep track of how many proposals were needed to get working numbers, highest count yet is <50
    count=0
    repeat {
      count=count+1
      # if(nrow(sp[1]>1 )){
      sp_prop <- lapply(1:length(sp), function(x) res[[x]][inf[[y]],])
      # }else{
      #   sp_prop <- lapply(1:length(sp), function(x) sp[[x]][y] )
      # }
      ##using sliding window proposal
      prop <-lapply(1:length(sp_prop), function(x) 10 * (stats::runif(1) - 0.5) + as.vector(t(sp_prop[[x]][1])) )
      
      for (x in 1:2){
        sp_prop[[x]][1] <- prop[[x]]
      }
      
      betas <- lapply(1:length(sp_prop), function(x) traits2coefs_sp( sp_prop[[x]] ) ) # Convert to beta coefficients
      
      yy = betas[[1]][,1] + betas[[1]][,2]*pa_data_inf[[y]]$X1+ betas[[1]][,3]*(pa_data_inf[[y]]$X1^2) + (betas[[2]][,1] + betas[[2]][,2]*pa_data_inf[[y]]$X2+ betas[[2]][,3]*(pa_data_inf[[y]]$X2^2) )
      
      presProb <- 1/(1+exp(-1*yy))
      
      likelihood= dbinom(pa_data[[y]]$y,1,prob=presProb,log=T)
      sumlikelihood= sum(likelihood)
      
      if (sumlikelihood!=-Inf){
        #paste("count",count)
        for (i in 1:length(sp)) {
          res_new[[i]][inf[[y]],]=sp_prop[[i]]
        }
        sum[[y]]<-sumlikelihood
        print(sum(unlist(sum)))
        
        ##brak only when the species reaches a non -inf likelihood
        print(count)
        break
      }
    }
  }
  
  betas <- lapply(1:length(res_new), function(x) traits2coefs(res_new[[x]])) # Convert to beta coefficients
  
  yy = lapply( 1:nrow(res_new[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1+ betas[[1]][i,3]*(pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2+ betas[[2]][i,3]*(pa_data[[i]]$X2^2) )))
  
  presProb <-lapply(1:nrow(res_new[[1]]), function(i) 1/(1+exp(-1*yy[[i]])) )
  
  yyp_mat<-lapply(1:nrow(res_new[[1]]), function(i) matrix(presProb[[i]], nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=T) )
  
  
  likelihood= lapply(1:nrow(res_new[[1]]), function(z) return(dbinom(pa_data[[z]]$y,1,prob=presProb[[z]],log=T)))
  sumlikelihood_new=lapply(1:length(likelihood), function(x)  sum(likelihood[[x]]))
  
  p_mat<-lapply(1:nrow(res_new[[1]]), function(i) matrix(likelihood[[i]], nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=T) )
  
  yyp_pred=yyp_mat
  lik_pred=p_mat
  
  if(plot==T){
    
    for ( i in 1:nrow(res[[1]]) )  {
      filled.contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), lik_pred[[i]], plot.title = {
        Arrows(res_new[[1]][[i,1]] +res_new[[1]][[i,2]],  res_new[[2]][[i,1]], res_new[[1]][[i,1]] - res_new[[1]][[i,2]],  res_new[[2]][[i,1]], code = 3, arr.length =.1  ,col=1, lty=1, lwd=2);
        Arrows(res_new[[1]][[i,1]], res_new[[2]][[i,1]] - res_new[[2]][[i,2]], res_new[[1]][[i,1]],  res_new[[2]][[i,1]] + res_new[[2]][[i,2]], code = 3, arr.length =.1  ,col=1, lty=1, lwd=2);
        contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred[[i]], col=1,  add=T);
        title(main = paste("likelihood surface fixed pres", "rep", "species", i))
      })
      
    }
    
    
  }
  res_ft_new<-lapply(1:length(res_new), function(x) t(apply(res_new[[x]], 1, forwardTransform1)))
  #lapply(1:length(res), function (x) rownames(res[[x]])<-tree$tip.label)
  #lapply(1:length(res), function (x) rownames(res_ft[[x]])<-tree$tip.label)
  
  res_td_new<-lapply(1:length(res_new), function(x) make.treedata(tree, res_ft_new[[x]])$dat)
  res_td_bt_new<-lapply(1:length(res_new), function(x) make.treedata(tree, res_new[[x]])$dat)
  
  return(list(start_pars_bt=res_new, start_pars_ft=res_ft_new, start_td=res_td_new, start_td_bt=res_td_bt_new, likelihoods=sumlikelihood_new))
}

MLglmStartpars<-function(species_data, tree, height=NULL, buffer=F){
  
  yy <- lapply(1:length(tree$tip.label), function(species) species_data[[species]]$y)
  X1 <- lapply(1:length(tree$tip.label), function(species) species_data[[species]]$X1 )
  X2 <- lapply(1:length(tree$tip.label), function(species) species_data[[species]]$X2 )
  starting.X1 = lapply(1:length(tree$tip.label), function(species) glm(yy[[species]] ~ X1[[species]] + I(X1[[species]]^2), family=binomial) ) #This is a null model
  
  
  #plot(X1[[62]], starting.X1[[62]]$fitted.values)
  
  #res_new[[2]][42,]
  
  starting.X2 = lapply(1:length(tree$tip.label), function(species) glm(yy[[species]] ~ X2[[species]] + I(X2[[species]]^2), family=binomial) )#This is a null model
  
  res_new<-list(
    do.call(rbind, lapply(1:length(tree$tip.label), function(species) coef2traits(starting.X1[[species]]$coefficients)))
    ,do.call(rbind, lapply(1:length(tree$tip.label), function(species) coef2traits(starting.X2[[species]]$coefficients)))
  )
  
  
  if(is.null(height)==F){
    res_new[[1]][,3]<-height
    res_new[[2]][,3]<-height
  }
  
  for (x in (1:length(res_new))){  rownames(res_new[[x]])<-tree$tip.label}
  
  if(buffer==T){
    
    for(sp in 1:length(tree$tip.label)){
      for ( trait in 1:2){
        
        if ( abs(res_new[[trait]][[sp,1]]) > 4 ){
          res_new[[trait]][[sp,1]]<-median(species_data[[sp]][[2+trait]][species_data[[sp]]$y==1])
        }
        
        if (abs(res_new[[trait]][[sp,2]])>4){
          res_new[[trait]][[sp,2]]<- diff(range(species_data[[sp]][[2+trait]][species_data[[sp]]$y==1]))
        }
      }
    }
    
  }
  
  res_ft_new<-lapply(1:length(res_new), function(x) t(apply(res_new[[x]], 1, forwardTransform1)))
  
  
  res_td_new<-lapply(1:length(res_new), function(x) make.treedata(tree, res_ft_new[[x]])$dat)
  res_td_bt_new<-lapply(1:length(res_new), function(x) make.treedata(tree, res_new[[x]])$dat)
  
  return(list(start_pars_bt=res_new, start_pars_ft=res_ft_new, start_td=res_td_new, start_td_bt=res_td_bt_new))
}

dataStartpars<-function(species_data, tree, height, randomize=F, rand_d=1){
  
  d=rand_d
  
  res_new<-list(
    do.call(rbind, lapply(1:length(tree$tip.label), function(species) c(1,1,1)))
    ,do.call(rbind, lapply(1:length(tree$tip.label), function(species) c(1,1,1)))
  )
  
  
  res_new[[1]][,3]<-height
  res_new[[2]][,3]<-height
  
  for (x in (1:length(res_new))){  rownames(res_new[[x]])<-tree$tip.label}
  
  if(randomize==F){
    
    for(sp in 1:length(tree$tip.label)){
      for ( trait in 1:2){
        
        res_new[[trait]][[sp,1]]<-median(species_data[[sp]][[2+trait]][species_data[[sp]]$y==1])
        
        
        res_new[[trait]][[sp,2]]<- diff(quantile(species_data[[sp]][[2+trait]][species_data[[sp]]$y==1], prob=c(.025,.975)))
        
        if(res_new[[trait]][[sp,2]]==0){
          res_new[[trait]][[sp,2]]=.05
        }
        
      }
    }
  }else{
    for(sp in 1:length(tree$tip.label)){
      for ( trait in 1:2){
        
        res_new[[trait]][[sp,1]]<- d * (stats::runif(1) - 0.5) + median(species_data[[sp]][[2+trait]][species_data[[sp]]$y==1])
        
        
        res_new[[trait]][[sp,2]]<-  d * (stats::runif(1) - 0.5) + diff(quantile(species_data[[sp]][[2+trait]][species_data[[sp]]$y==1], prob=c(.025,.975)))
      }
    }
  }
  
  
  
  
  
  res_ft_new <-lapply(1:length(res_new), function(x) t(apply(res_new[[x]], 1, forwardTransform1)))
  
  
  res_td_new <-lapply(1:length(res_new), function(x) make.treedata(tree, res_ft_new[[x]])$dat)
  res_td_bt_new <-lapply(1:length(res_new), function(x) make.treedata(tree, res_new[[x]])$dat)
  
  return(list(start_pars_bt=res_new, start_pars_ft=res_ft_new, start_td=res_td_new, start_td_bt=res_td_bt_new))
}

find_good_start_pars=function(Prior_object, tree, dist="norm",  ntries=1000){
  rep=0
  
  glm_heights=lapply(1:length(Prior_object), function(pred) Prior_object[[pred]]$pars$heights)
  
  repeat{
    print(rep)
    rep=rep+1
    #print
    startPars_scaled <- priorSim_pars(Prior_object, tree, dist=Prior_object[[1]]$pars$den.mu,hard_coded_heights = glm_heights )
    #startPars_scaled$sim_da$sim_dat_ft=GLM_only_ml$start_pars_ft
    #startPars_scaled$sim_da$sim_dat_bt=GLM_only_ml$start_pars_bt
    #startPars_scaled$sim_da$sim_td=GLM_only_ml$start_td
    #startPars_scaled$sim_da$sim_td_bt=GLM_only_ml$start_td_bt
    
    
    
    #break when no tip has a -inf likelihood
    
    if(length(findBadStart(res=startPars_scaled$sim_dat$sim_dat_bt, pa_data=data_final, plot=T))==0){
      break
      return(startPars_scaled)
    }
    
    if(rep>ntries){
      break
      return("could not find starting values in ntries")
    }
    
    
  }
  
}


MLglmStartpars_general = function(species_data, tree, height = NULL, buffer = FALSE) {
  n_species <- length(tree$tip.label)
  
  # Detect how many X variables are present (assuming structure of species_data[[i]]$X1, $X2, ..., $Xn)
  x_names <- names(species_data[[1]])
  x_vars <- grep("^X[0-9]+$", x_names, value = TRUE)
  nx <- length(x_vars)
  
  # Extract y and each X for all species
  yy <- lapply(1:n_species, function(i) species_data[[i]]$y)
  XX <- lapply(x_vars, function(xn) lapply(1:n_species, function(i) species_data[[i]][[xn]]))
  
  # Fit GLMs per predictor per species
  starting_models <- lapply(1:nx, function(xi) {
    lapply(1:n_species, function(i) {
      glm(yy[[i]] ~ XX[[xi]][[i]] + I(XX[[xi]][[i]]^2), family = binomial)
    })
  })
  
  # Extract and convert coefficients
  res_new <- lapply(1:nx, function(xi) {
    do.call(rbind, lapply(1:n_species, function(i) coef2traits(starting_models[[xi]][[i]]$coefficients)))
  })
  
  # Add height if specified
  if (!is.null(height)) {
    for (x in 1:nx) {
      res_new[[x]][, 3] <- height
    }
  }
  
  # Add rownames to result
  for (x in 1:nx) {
    rownames(res_new[[x]]) <- tree$tip.label
  }
  
  # Apply buffer logic if enabled
  if (buffer) {
    for (i in 1:n_species) {
      for (xi in 1:nx) {
        # Assume X variables are stored after y in the list, i.e., y, X1, X2, ..., Xn
        X_vals <- species_data[[i]][[x_vars[xi]]][species_data[[i]]$y == 1]
        if (abs(res_new[[xi]][i, 1]) > 4) {
          res_new[[xi]][i, 1] <- median(X_vals)
        }
        if (abs(res_new[[xi]][i, 2]) > 4) {
          res_new[[xi]][i, 2] <- diff(range(X_vals))
        }
      }
    }
  }
  
  # Apply forward transform and convert to treedata
  res_ft_new <- lapply(res_new, function(mat) t(apply(mat, 1, forwardTransform1)))
  res_td_new <- lapply(res_ft_new, function(mat) make.treedata(tree, mat)$dat)
  res_td_bt_new <- lapply(res_new, function(mat) make.treedata(tree, mat)$dat)
  
  return(list(
    start_pars_bt = res_new,
    start_pars_ft = res_ft_new,
    start_td = res_td_new,
    start_td_bt = res_td_bt_new
  ))
}


get_starting_values = function(Prior_scale, reps_before_POE=1000){
  
  GLM_only_ml<-suppressWarnings(MLglmStartpars_general(species_data = data_obj$data,tree = tree, height = NULL))
  
  heights_glm<- lapply(GLM_only_ml$start_pars_bt, function(pred) pred[,3])
  
  #To get start values we simulate parameters of our model, including individual species response curves and mvBM pars, however because we only can set what the mvBM pars are and the tip responses are subsequently simulated, we dont have a lot of control over what the tip values will be and some potential starting values could be so far off they
  
  #```{r}
  # Simulate mvBM and starting individual species responses
  startPars_scaled=list()
  rep=0
  #startPars_scaled=GLM_only_ml$
  
  bad_names=list()
  
  reps=1
  
  ##repeatedly simulate response curve characters for each species until you get a starting set that gives you a tractable loglikleihood
  # reps indiciates number of independent MCMC's
  for (i in 1:reps){
    Prior_scale
    rep=0
    
    repeat{
      print(rep)
      rep=rep+1
      #print
      
      #for the firsts 10000 tries just resmaple ALL response curves
      if(rep<5){
        startPars_scaled[[i]] <- priorSim_pars(Prior_scale, tree, dist="norm",hard_coded_heights = heights_glm)
      }else {
        
        #after 1000 iterations, sample only for the curves that are still NA's and repeast until none of them are NA's
        startPars_scaled_draw <- priorSim_pars(Prior_scale, tree, dist="norm",hard_coded_heights = heights_glm)
        #startPars_scaled[[i]]$sim_dat$sim_dat_bt[[1]][bad_idx,1:2]=startPars_scaled_draw$sim_dat$sim_dat_bt[[1]][bad_idx,1:2]
        #startPars_scaled[[i]]$sim_dat$sim_dat_bt[[2]][bad_idx,1:2]=startPars_scaled_draw$sim_dat$sim_dat_bt[[2]][bad_idx,1:2]
        #
        #startPars_scaled[[i]]$sim_dat$sim_dat_bt[[1]][bad_idx,1:2]=startPars_scaled_draw$sim_dat$sim_dat_bt[[1]][bad_idx,1:2]
        #startPars_scaled[[i]]$sim_dat$sim_dat_bt[[2]][bad_idx,1:2]=startPars_scaled_draw$sim_dat$sim_dat_bt[[2]][bad_idx,1:2]
        #
        #startPars_scaled[[i]]$sim_dat$sim_td_bt[[1]][,1] = startPars_scaled[[i]]$sim_dat$sim_dat_bt[[1]][,1]
        #startPars_scaled[[i]]$sim_dat$sim_td_bt[[1]][,2] = startPars_scaled[[i]]$sim_dat$sim_dat_bt[[1]][,2]
        #
        #
        #startPars_scaled[[i]]$sim_dat$sim_td_bt[[2]][,1] = startPars_scaled[[i]]$sim_dat$sim_dat_bt[[2]][,1]
        #startPars_scaled[[i]]$sim_dat$sim_td_bt[[2]][,2] = startPars_scaled[[i]]$sim_dat$sim_dat_bt[[2]][,2]
        if (rep < reps_before_POE) {
          startPars_scaled[[i]] <- priorSim_pars(Prior_scale, tree, dist = "norm", hard_coded_heights = heights_glm)
        } else {
          startPars_scaled_draw <- priorSim_pars(Prior_scale, tree, dist = "norm", hard_coded_heights = heights_glm)
          
          for (k in seq_along(startPars_scaled[[i]]$sim_dat$sim_dat_bt)) {
            startPars_scaled[[i]]$sim_dat$sim_dat_bt[[k]][bad_idx, 1:2] <- startPars_scaled_draw$sim_dat$sim_dat_bt[[k]][bad_idx, 1:2]
            startPars_scaled[[i]]$sim_dat$sim_td_bt[[k]][, 1:2] <- startPars_scaled[[i]]$sim_dat$sim_dat_bt[[k]][, 1:2]
          }
        }
        
        
      }
      #startPars_scaled$sim_da$sim_dat_bt=GLM_only_ml$start_pars_bt
      #startPars_scaled$sim_da$sim_td=GLM_only_ml$start_td
      #startPars_scaled$sim_da$sim_td_bt=GLM_only_ml$start_td_bt
      #break when no tip has a -inf likelihood
      
      #book keep which species are stille bad curves (both in terms of species name and row index in the response curve dataframe)
      bad_names[[rep]]=names(heights_glm[[1]])[findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data=data_final, plot=T)]
      bad_idx=findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data=data_final, plot=T)
      
      #if none are bad starting values, break the loop, you are good to go! if not, keep going!
      if(length(findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data=data_final, plot=T))==0){break}
      
    }
    #
    findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data=data_final, plot=F)
    
    
    #startPars_scaled[[i]]$sim_dat$sim_dat_bt[[1]][,2] ==startPars_scaled[[i]]$sim_dat$sim_td_bt[[1]][,2]
  }
  
  return(startPars_scaled)
}


sim_df2listoflists_fn<-  function(sim_dat, tree, Multi_Pred=T, k){####
  #produces a list of lists of traits for each species
  #part of sim_data_fn
  #provide simulated dfata from MvMorph::mvSIM
  #must use  untransformed values
  
  Multi_Pred=T
  
  
  if(Multi_Pred==T){
    dat= sim_dat
    dat_final<-list()
    df_final<-list()
    
    for (i in 1:length(dat)){
      
      dat0=t(apply(dat[[i]], 1, backTransform1))
      
      dat1=as.data.frame(dat0)
      
      dat1=as.data.frame(dat1)
      
      dat_final[[i]]<- tibble::rownames_to_column(dat1, "VALUE")
      
      #putting names to list so simulation code is legible, df is the final object
      df_list <- lapply(1:nrow(dat_final[[i]]), function(j) as.list(dat_final[[i]][j,]))
      df_sp<- lapply(df_list, function(x) {
        names(x)[grep("VALUE", names(x))] <- "species"
        x})
      df_theta<- lapply(df_sp, function(x) {
        names(x)[grep("V1", names(x))] <- "theta"
        x})
      df_width<- lapply(df_theta, function(x) {
        names(x)[grep("V2", names(x))] <- "width"
        x})
      df<- lapply(df_width, function(x) {
        names(x)[grep("V3", names(x))] <- "height"
        x})
      #simulating presence/absence data
      df_final[[i]]<-df
      
      
      # View(df_final)
    }
    # return(df_final)
    
  } else {
    
    
    #simulate data for 3 traits given R and Theta
    dat0=t(apply(dat, 1, backTransform1))
    
    dat1=as.data.frame(dat0)
    
    dat1=as.data.frame(dat1)
    dat2<- tibble::rownames_to_column(dat1, "VALUE")
    
    #putting names to list so simulation code is legible, df is the final object
    df_list <- lapply(1:nrow(dat2), function(j) as.list(dat2[j,]))
    df_sp<- lapply(df_list, function(x) {
      names(x)[grep("VALUE", names(x))] <- "species"
      x})
    df_theta<- lapply(df_sp, function(x) {
      names(x)[grep("V1", names(x))] <- "theta"
      x})
    df_width<- lapply(df_theta, function(x) {
      names(x)[grep("V2", names(x))] <- "width"
      x})
    df<- lapply(df_width, function(x) {
      names(x)[grep("V3", names(x))] <- "height"
      x})
    #simulating presence/absence data
  }
  #}
  #View(df)
  #######generate betas{#####
  #gens regression coefficients and sim temp data for existing list of lists
  #part of sim_data_fn
  #needs output from sim_traits_fn (list of lists of trait data for each species)
  if (Multi_Pred==T) {
    
    df<-df_final
    
    for (x in (1:length(df))){
      #x=1
      
      
      for (i in (1:length(df[[x]]))){
        #i=1
        #real_H
        df[[x]][[i]]$H= log((1/(df[[x]][[i]]$height))-1)
        
        df[[x]][[i]]$b0 <- (((df[[x]][[i]]$theta)^2)/((df[[x]][[i]]$width)^2))*k[[x]] + (df[[x]][[i]]$H)*((((df[[x]][[i]]$theta)^2)/((df[[x]][[i]]$width)^2))-1)
        df[[x]][[i]]$b1 <- (-2*((df[[x]][[i]]$H)+k[[x]])*(df[[x]][[i]]$theta))/(df[[x]][[i]]$width)^2
        df[[x]][[i]]$b2<- (1/(df[[x]][[i]]$width)^2)*((df[[x]][[i]]$H)+k[[x]])
        
      }
      
    }
    #View(df)
    return(df)
  } else{
    for (i in (1:length(df))){
      
      #real_H
      df[[i]]$H= log((1/(df[[i]]$height))-1)
      
      df[[i]]$b0 <- (((df[[i]]$theta)^2)/((df[[i]]$width)^2))*k + (df[[i]]$H)*((((df[[i]]$theta)^2)/((df[[i]]$width)^2))-1)
      df[[i]]$b1 <- (-2*((df[[i]]$H)+k)*(df[[i]]$theta))/(df[[i]]$width)^2
      df[[i]]$b2<- (1/(df[[i]]$width)^2)*((df[[i]]$H)+k)
    }
    return(df)
  }
}


simPA<-function(res, tree, span, grid_size, simMiss=F, nMiss=1){
  
  
  if(simMiss==T){
    miss<-sample(1:length(tree$tip.label), nMiss)
  }
  
  
  #generate bound for PA/ sim for each species
  max=lapply(1:length(res), function(pred) lapply( 1:nrow(res[[1]]), function(sp) round(res[[pred]][[sp,1]]+(res[[pred]][[sp,2]]), 3)))
  min=lapply(1:length(res), function(pred) lapply( 1:nrow(res[[1]]), function(sp) round(res[[1]][[sp,1]]-(res[[1]][[sp,2]]), 3)))
  
  #generate predictor vectors using seq based on what the expanded.grid output length will be (grid size being nth rooted based on # of pred)
  Pred_full<- lapply( 1:nrow(res[[1]]), function(sp) lapply(1:length(res), function(pred) seq(min[[pred]][[sp]]-span,max[[pred]][[sp]]+span, length.out = (grid_size^ (1 / length(res))) )))
  #Y_full<-lapply( 1:nrow(res[[1]]), function(sp) seq(ymin[[sp]]-span,ymax[[sp]]+span, length.out = sqrt(grid_size)))
  
  #generate predictor grid for each sp
  full_grid=lapply( 1:nrow(res[[1]]), function(sp) as.matrix(expand.grid( Pred_full[[sp]] )) )
  
  betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
  
  #yy_prelogit  = lapply( 1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*full_grid[[i]][,1] + betas[[1]][i,3]*(full_grid[[i]][,1] ^2) + (betas[[2]][i,1] + betas[[2]][i,2]*full_grid[[i]][,2]+ betas[[2]][i,3]*(full_grid[[i]][,2] ^2) )) )
  
  #made the betas iterative so that we can flexibly simulate multiple pred
  
  yy_prelogit_sep  = lapply( 1:nrow(res[[1]]), function(sp)  lapply(1:length(res), function(pred) (betas[[pred]][sp,1] + betas[[pred]][sp,2]*full_grid[[sp]][,pred]+ betas[[pred]][sp,3]*(full_grid[[sp]][,pred]^2)  ) ) )
  
  yy_prelogit= lapply(1:length(yy_prelogit_sep), function(sp) rowSums(matrix(unlist(yy_prelogit_sep[[sp]]), ncol=length(yy_prelogit_sep[[sp]]), byrow=F)))
  
  presProb <-lapply(1:nrow(res[[1]]), function(i) 1/(1+exp(-1*yy_prelogit[[i]])) ) #convert to probability using logit link
  
  y <- lapply(1:nrow(res[[1]]), function(i) rbinom(n=length(presProb[[i]]), size=1, prob=presProb[[i]]) )
  
  names(y)<- lapply(1:nrow(res[[1]]), function(i)  "y")
  
  species_data=list()
  
  for (sp in  1:nrow(res[[1]])){
    
    species_data[[sp]]<-list(
      species= tree$tip.label[[sp]],
      y= y[[sp]]
    )
    
    for(pred in 1:ncol(full_grid[[sp]])){
      
      species_data[[sp]][[paste("X", pred, sep="") ]] <-full_grid[[sp]][,pred]
      #View(species_data)
      
      #species_data[[sp]][[length(species_data[[sp]])+1]] <- list((paste("X", pred, sep="") )=full_grid[[sp]][,pred])
      
      #names(species_data[[sp]][length(species_data[[sp]])])<- paste("X", pred, sep="")
    }
  }
  
  
  if(simMiss==T){
    
    for (sp in miss){
      
      species_data[[sp]] <- list(species= tree$tip.label[[sp]],
                                 y= NA
      )
      for(pred in 1:ncol(full_grid[[sp]])){
        
        species_data[[sp]][[paste("X", pred, sep="") ]] <- NA
      }
      
    }
    return(list(species_data=species_data, missing_sp= sort(miss)))
    
  }
  
  
  
  return(list(species_data=species_data))
  
}
