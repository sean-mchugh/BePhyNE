#require(readr)
#require(tibble)
#require(MCMCpack)
#require(devtools)
#require(coda)
#require(devtools)
#require(coda)
#require(ape)
#require(truncdist)
#require(geiger)
#require(phytools)
#require(tidyr)
#require(ratematrix)
#require(mvMORPH)
#require(dplyr)
#require(mvtnorm)
#require(MultiRNG)
#require(Rphylopars)
#require(Rcpp)
#require(doParallel)
#require(corpcor)
#require(Matrix)
#require(treeplyr)
##require(bayou)
#require(crayon)
#require(shape)
#require(scales)
#require(phytools)
#require(robustbase)



{

  ###Utils##################################################################################################################################################################

  data_df2list=function(df, tree){

    #```{r}
    #
    sp_col=(1:ncol(df))[names(df)=="species"]

    data_final<- lapply(split(df[sp_col:ncol(df)],df$species), as.list)
    #replace vector of species names with a single name
    for (sp in 1:length(data_final)) {

      data_final[[sp]]$species= data_final[[sp]]$species[[1]]

    }

    # match data_final species order to tree tip order
    data_final=data_final[order(match(unlist(lapply(data_final, function(sp) sp$species)), tree$tip.label))]

    #```

    return(data_final)
  }

  make.treedata =function (tree, data, name_column = "detect", as.is = FALSE)
  {
    #from Josef Uyeda's treeplyer package loaded here to avoid having to install treeplyer
    if (class(tree) != "phylo")
      stop("tree must be of class 'phylo'")
    if (is.vector(data)) {
      data <- as.matrix(data)
      colnames(data) <- "trait"
    }
    if (is.null(colnames(data))) {
      colnames(data) <- paste("trait", 1:ncol(data), sep = "")
    }
    coln <- colnames(data)
    if (name_column == "detect") {
      if (is.null(rownames(data))) {
        tmp.df <- data.frame(data)
        offset <- 0
      }
      else {
        tmp.df <- data.frame(rownames(data), data)
        offset <- 1
      }
      matches <- sapply(tmp.df, function(x) sum(x %in% tree$tip.label))
      if (all(matches == 0))
        stop("No matching names found between data and tree")
      name_column <- which(matches == max(matches)) - offset
    }
    else {
      if (is.character(name_column)) {
        name_column <- which(name_column == coln)[1]
      }
    }
    dat <- as_tibble(as.data.frame(lapply(1:ncol(data), function(x) type.convert(apply(data[,
                                                                                            x, drop = FALSE], 1, as.character), as.is = as.is))))
    colnames(dat) <- coln
    if (name_column == 0) {
      clnm <- colnames(dat)
      dat <- dat[, clnm, drop = FALSE]
      dat.label <- as.character(rownames(data))
    }
    else {
      if (is.numeric(name_column)) {
        clnm <- (1:ncol(data))[-name_column]
      }
      else {
        clnm <- colnames(dat)[-which(colnames(dat) == name_column)]
      }
      dat <- dat[, clnm, drop = FALSE]
      dat.label <- as.character(as.data.frame(data)[[name_column]])
    }
    data_not_tree <- setdiff(dat.label, tree$tip.label)
    tree_not_data <- setdiff(tree$tip.label, dat.label)
    phy <- drop.tip(tree, tree_not_data)
    dat <- filter(dat, dat.label %in% phy$tip.label)
    dat.label <- dat.label[dat.label %in% phy$tip.label]
    if (any(duplicated(dat.label))) {
      warning("Duplicated data in dataset, selecting first unique entry for each species")
      dat <- filter(dat, !duplicated(dat.label))
      dat.label <- dat.label[!duplicated(dat.label)]
    }
    ...my.order... <- match(dat.label, phy$tip.label)
    dat <- arrange(dat, ...my.order...)
    td <- list(phy = phy, dat = dat)
    class(td) <- c("treedata", "list")
    attributes(td)$tip.label <- phy$tip.label
    attributes(td)$dropped <- list(dropped_from_tree = data_not_tree,
                                   dropped_from_data = tree_not_data)
    return(td)
  }


  makeTransparent<- function (someColor, alpha = 100)
  {
    newColor <- col2rgb(someColor)
    apply(newColor, 2, function(curcoldata) {
      rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3],
          alpha = alpha, maxColorValue = 255)
    })
  }

  ###transformations##################################################################################################################################################################

  traits2coefs <- function(traits, v=0.05){


    if (ncol(as.matrix(traits))==1){
      theta=traits[[1]]
      w=traits[[2]]
      height=traits[[3]]
    } else{
      theta=traits[,1]
      w=traits[,2]
      height=traits[,3]

    }

    #theta=traits[,1]
    #w=traits[,2]
    #height=traits[,3]

    H=log((1/height)-1)
    k=log(v/(1-v))
    b0<- ((theta^2)/(w^2))*k + H*(((theta^2)/(w^2))-1)
    b1<- (-2*(H+k)*theta)/w^2
    b2<- (1/w^2)*(H+k)
    return(data.frame(b0=b0,b1=b1, b2=b2))
  }

  traits2coefs_sp<-function(traits, v=0.05){
    theta=traits[1]
    w=traits[2]
    height=traits[3]
    H=log((1/height)-1)
    k=log(v/(1-v))
    b0<- ((theta^2)/(w^2))*k + H*(((theta^2)/(w^2))-1)
    b1<- (-2*(H+k)*theta)/w^2
    b2<- (1/w^2)*(H+k)
    return(data.frame(b0=b0,b1=b1, b2=b2))
  }


  coef2traits<-function(coefs, v=0.05){
    B0=coefs[[1]]
    B1=coefs[[2]]
    B2=coefs[[3]]

    theta=B1*-1/(2*B2)
    height=1 / (1 + exp( B1^2 / (4*B2) - B0) )
    width = -sqrt( B1^2 - 4*B2 * (B0 - log( (v)/(1-v) ) ) ) / (2*B2)

    return(c(theta, width, height))

  }


  forwardTransform1 <- function(x){
    newx <- x
    newx[2] <- log(x[2])
    Y=((x[3]-.05)/.95)
    newx[3] <- -1*log(Y/(1-Y))
    return(newx)
  }


  backTransform1 <- function(x){
    newx <- x
    newx[2] <- exp(x[2])
    newx[3] <- 0.05 + 0.95*(exp(-1*x[3])/(1+exp(-1*x[3])))
    return(newx)

  }


  unscale<- function(dat, scale_atr) {
    unscaled<- lapply(1:length(median), function(pred) cbind(dat[[pred]][,1]*scale_atr$scale[[pred]]+scale_atr$center[[pred]],
                                                             dat[[pred]][,2]*scale_atr$scale[[pred]],
                                                             dat[[pred]][,3])
    )

    return(unscaled)


  }

  C2K<-function(C){

    K = C  + 273.15

    return(K)
  }


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


  #######Prior#############################################################################################################################################################
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



  makePrior_ENE<- function(r, p, den.mu="unif", par.mu, den.sd="unif", den.heights="unif", par.sd, heights=NULL, unif.corr=TRUE, Sigma=NULL, nu=NULL, plot=T){


    par.mu.c <- par.mu[1,]
    par.mu.w <- par.mu[2,]

    if(plot==T){
      par(mfrow = c(2, 2))
    }


    ## Make a warning if 'unif.corr' is TRUE and 'Sigma' or 'nu' has values.
    if( unif.corr ){
      if( !is.null(Sigma) | !is.null(nu) ){
        warning("Found values for 'Sigma' and 'nu', but 'unif.corr=TRUE'. Using uniform correlations.")
      }
    }

    ## Save all parameters to use in subsequent functions.
    pars <- list()
    pars$r <- r ## Number of traits
    pars$p <- p ## Number of rate regimes.
    pars$den.mu <- den.mu
    pars$par.mu <-par.mu
    pars$par.mu.c <- par.mu[1,]
    pars$par.mu.w <- par.mu[2,]
    pars$den.sd <- den.sd
    pars$par.sd <- par.sd
    pars$unif.corr <- unif.corr
    pars$Sigma <- Sigma
    pars$nu <- nu
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
      mn<- function(x) sum( c(stats::dnorm(x[1], mean=par.mu.c[[1]], sd=par.mu.c[[2]], log=TRUE),
                              stats::dlnorm(x[2],mean=par.mu.w[[1]], sd=par.mu.w[[2]], log=TRUE)
      ) )

      if(plot==T){
        plot(density(stats::rnorm(10000, mean=par.mu[1,1], sd=par.mu[1,2]) ), main=paste("A_c") )
        plot(density(stats::rlnorm(10000, meanlog=par.mu[2,1], sdlog=par.mu[2,2]) ), main=paste("A_w")  )
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

    if(is.null(heights)==T){

      height_priors=NULL

    }else{



      if(den.heights=="unif"){

        height_priors<-lapply(1:length(heights), function(sp){ force(heights[[sp]]);
          (function(x) dunif(x, min = -4.5, max = 5, log=T) )  })


      }else{


        height_priors<-list()



        Y=((heights-.05)/.95)
        FT_heights<- -1*log(Y/(1-Y))

        #height_priors=list()




        #making a list of functions rather than one big function so we can call the specific function we want only when the move on that par happens,
        #do we really need to calculate the prior for EVERY sp EVERY iteration...nah thats gonna slow the whole mcmc down so instead lets just do it when there is a proposal on that height
        height_priors<-lapply(1:length(FT_heights), function(sp){ force(FT_heights[[sp]]);
          (function(x) dnorm(x, mean=FT_heights[[sp]], 0.15, log=T) )  })

      }

    }


    res <- list(mean.prior=mn, corr.prior=corr, sd.prior=sd,  pars=pars, heights.prior=height_priors)
    class( res ) <- "prior_function"
    return( res )
  }

  makePrior_ENE<- function(r, p, den.mu="unif", par.mu, den.sd="unif", den.heights="unif", par.sd, heights=NULL, unif.corr=TRUE, Sigma=NULL, nu=NULL, plot=T){


    par.mu.c <- par.mu[1,]
    par.mu.w <- par.mu[2,]

    if(plot==T){
      par(mfrow = c(2, 2))
    }


    ## Make a warning if 'unif.corr' is TRUE and 'Sigma' or 'nu' has values.
    if( unif.corr ){
      if( !is.null(Sigma) | !is.null(nu) ){
        warning("Found values for 'Sigma' and 'nu', but 'unif.corr=TRUE'. Using uniform correlations.")
      }
    }

    ## Save all parameters to use in subsequent functions.
    pars <- list()
    pars$r <- r ## Number of traits
    pars$p <- p ## Number of rate regimes.
    pars$den.mu <- den.mu
    pars$par.mu <-par.mu
    pars$par.mu.c <- par.mu[1,]
    pars$par.mu.w <- par.mu[2,]
    pars$den.sd <- den.sd
    pars$par.sd <- par.sd
    pars$unif.corr <- unif.corr
    pars$Sigma <- Sigma
    pars$nu <- nu
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
      mn<- function(x) sum( c(stats::dnorm(x[1], mean=par.mu.c[[1]], sd=par.mu.c[[2]], log=TRUE),
                              stats::dlnorm(x[2],mean=par.mu.w[[1]], sd=par.mu.w[[2]], log=TRUE)
      ) )

      if(plot==T){
        plot(density(stats::rnorm(10000, mean=par.mu[1,1], sd=par.mu[1,2]) ), main=paste("A_c") )
        plot(density(stats::rlnorm(10000, meanlog=par.mu[2,1], sdlog=par.mu[2,2]) ), main=paste("A_w")  )
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

    if(is.null(heights)==T){

      height_priors=NULL

    }else{



      if(den.heights=="unif"){

        height_priors<-lapply(1:length(heights), function(sp){ force(heights[[sp]]);
          (function(x) dunif(x, min = -4.5, max = 5, log=T) )  })


      }else{


        height_priors<-list()



        Y=((heights-.05)/.95)
        FT_heights<- -1*log(Y/(1-Y))

        #height_priors=list()




        #making a list of functions rather than one big function so we can call the specific function we want only when the move on that par happens,
        #do we really need to calculate the prior for EVERY sp EVERY iteration...nah thats gonna slow the whole mcmc down so instead lets just do it when there is a proposal on that height
        height_priors<-lapply(1:length(FT_heights), function(sp){ force(FT_heights[[sp]]);
          (function(x) dnorm(x, mean=FT_heights[[sp]], 0.15, log=T) )  })

      }

    }














    res <- list(mean.prior=mn, corr.prior=corr, sd.prior=sd,  pars=pars, heights.prior=height_priors)
    class( res ) <- "prior_function"
    return( res )
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

        } else{

          mn<- function(x) sum( c(stats::dnorm(x[1], mean=par.mu.c[[1]], sd=par.mu.c[[2]], log=TRUE),
                                  stats::dlnorm(x[2],mean=par.mu.w[[1]], sd=par.mu.w[[2]], log=TRUE)
          ) )

          if(plot==T){
            plot(density(stats::rnorm(10000, mean=par.mu[1,1], sd=par.mu[1,2]) ), main=paste("A_c") )
            plot(density(stats::rlnorm(10000, meanlog=par.mu[2,1], sdlog=par.mu[2,2]) ), main=paste("A_w")  )

          }


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
          plot(density(stats::rlnorm(10000, meanlog=par.sd[3,1], sdlog=par.sd[3,2]) ), main=paste("r_sd_h")  )


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



      Y=((heights-.05)/.95)
      FT_heights<- -1*log(Y/(1-Y))

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


  #######proposal functions#########################################################################################################################################################################


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



  #likelihood_fns#####################################################################################################################################################################


  calc_pl <- function(b0, b1, b2, data, calc_pred=T){
    #calc probability of occurence from P/A data for one species
    pl=(b2*as.numeric(data)^2+b1*as.numeric(data) + b0)
    if (calc_pred==T){
      pred=exp(pl)/(1+exp(pl))
      return(pred)
    }
    else{
      return(pl)
    }
  }




  #fixed missing dat bug

  GLM_Likelihood_MV_fn <- function(res, full_pa_data, k){
    #td<-prop_traits_only

    #identify missing by if the first pa_point is NA (Do not have NAs unless there is no data for that sp)
    miss=unlist(lapply(1:length(full_pa_data), function(sp) if(is.na(full_pa_data[[sp]]$y[[1]])==T){sp} ))
    #save vector of which list entries contain data and which dont
    sp_dat<-(1:length(full_pa_data))[(1:length(full_pa_data))!=miss]





    if(is.null(miss)==T){
      sp_dat<-1:length(full_pa_data)
    }





    #check to see if the only sp is this missing data

    if(length(sp_dat)==0& length(miss)>0){
      sumlikelihood=list()
      sumlikelihood[[1]]<-0
      return(sumlikelihood)
    }

    #print(paste("miss",length(miss)))
    #print(paste("spdat",length(sp_dat)))

    if(is.null(miss)==F){
      #iterate through all species that dont have missing data
      pa_data <- lapply((1:length(full_pa_data))[(1:length(full_pa_data))!=miss], function(msp) full_pa_data[[msp]])
    }else{
      pa_data<-full_pa_data
    }

    #if proposing on only one tip, and that tips missing, you are going to get an error because pa_data will be empty

    #print(paste("length(pa_data)", length(pa_data)))
    #print(class(res[[1]]))
    X<-lapply(1:length(res), function(pred) lapply(1:length(full_pa_data), function(sp)   full_pa_data[[sp]][[2+pred]]) )

    # print(as.matrix(res[[1]]))
    betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients


    # yy_prelogit_sep  = lapply( 1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1+ betas[[1]][i,3]*(pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2+ betas[[2]][i,3]*(pa_data[[i]]$X2^2) )) )
    #with betas[[1]] error
    #yy_prelogit_sep  = lapply( 1:nrow(res[[1]]), function(sp)  lapply(1:length(res), function(pred) (betas[[1]][sp,1] + betas[[1]][sp,2]*X[[pred]][[sp]]+ betas[[1]][sp,3]*(X[[pred]][[sp]]^2)  ) ) )


    yy_prelogit_sep  = lapply(sp_dat, function(sp)  lapply(1:length(res), function(pred) (betas[[pred]][sp,1] + betas[[pred]][sp,2]*X[[pred]][[sp]]+ betas[[pred]][sp,3]*(X[[pred]][[sp]]^2)  ) ) )

    yy_prelogit= lapply(1:length(yy_prelogit_sep), function(sp) rowSums(matrix(unlist(yy_prelogit_sep[[sp]]), ncol=length(yy_prelogit_sep[[sp]]), byrow=F)))



    presProb <-lapply(1:length(yy_prelogit), function(i) 1/(1+exp(-1*yy_prelogit[[i]])) ) #convert to probability using logit link
    #likelihood=dnorm(y,mean=pred,sd=sd, log=T)
    likelihood= lapply(1:length(presProb), function(z) return(dbinom(pa_data[[z]]$y,1,prob=presProb[[z]],log=T)))
    sumlikelihood=lapply(1:length(likelihood), function(x)  sum(likelihood[[x]]))

    if(is.null(miss)==F){
      sp_dat_lnl<-list()
      for (sp in 1:length(sp_dat)){
        sp_dat_lnl[[sp_dat[[sp]]]]<-sumlikelihood[[sp]]

      }
      for (sp in miss){
        sp_dat_lnl[[sp]]<-0

      }
      return(sp_dat_lnl)
    }



    return(sumlikelihood)
  }



  lnL_ratematrix <- function (proposal, tree, root, R, H_fixed=T) {
    if (H_fixed==T){
      X.c <- sapply(1:(ncol(proposal)-1), function(x) proposal[, x] - root[x])

    } else{
      X.c <- sapply(1:ncol(proposal), function(x) proposal[, x] - root[x])
    }
    P <- data.frame(X.c, row.names = tree$tip.label)
    comp <- phylolm::three.point.compute(tree, P = P, Q = NULL)
    xiVx <- sum(comp$PP * chol2inv(chol(R)))
    detV <- (ncol(proposal)-1) * comp$logd + nrow(proposal) * determinant(R)$modulus[1]
    logl <- -0.5 * ((ncol(proposal)-1) * nrow(proposal) * log(2 * pi) + detV + xiVx)
    return(logl)

  }

  lnL_ratematrix <- function (proposal, tree, root, R) {
    X.c <- sapply(1:length(root), function(x) proposal[, x] - root[x])

    P <- data.frame(X.c, row.names = tree$tip.label)
    comp <- phylolm::three.point.compute(tree, P = P, Q = NULL)
    xiVx <- sum(comp$PP * chol2inv(chol(R)))
    detV <- (length(root)) * comp$logd + nrow(proposal) * determinant(R)$modulus[1]
    logl <- -0.5 * ((length(root)) * nrow(proposal) * log(2 * pi) + detV + xiVx)
    return(logl)

  }



  hastingsDensity <- function(curr.vcv, prop.vcv, p, v){
    ## Calculates the log Hastings ratio by dividing p(curr|prop)/p(prop|curr) densities.
    center.curr <- (v-p-1) * curr.vcv
    center.prop <- (v-p-1) * prop.vcv
    hh <- logDensityIWish(curr.vcv, v, center.prop) - logDensityIWish(prop.vcv, v, center.curr)
    return(hh)
  }


  ###calculating total likelihood #########################################################################################################################################################################
  acceptance_ratio_MV_fun<-function(proposal, current_vals, pa_data, tree, prior, prop, k, prior_only=F, glm_only=F, H_fixed=T){

    dat.cur=lapply(1:length(current_vals[[1]][[1]]), function (x) current_vals[[1]][[1]][[x]]$dat)
    dat.cur.bt=lapply(1:length(dat.cur), function(x) t(apply(dat.cur[[x]], 1, backTransform1)))

    dat.prop=lapply(1:length(proposal), function (x) proposal[[x]]$td$dat)
    dat.prop.bt=lapply(1:length(dat.prop), function(x) t(apply(dat.prop[[x]], 1, backTransform1)))

    A.prop=lapply(1:length(proposal), function (x) proposal[[x]]$A)
    A.prop.bt=lapply(1:length(A.prop), function(x) backTransform1(A.prop[[x]])[1:length(A.prop[[x]])])

    A.cur=lapply(1:length(current_vals[[1]][[1]]), function (x) current_vals[[1]][[1]][[x]]$A)
    A.cur.bt=lapply(1:length(A.cur), function(x) backTransform1(A.cur[[x]])[1:length(A.cur[[x]])])

    R_cor.prop=lapply(1:length(proposal), function (x) proposal[[x]]$R_cor)
    R_sd.prop=lapply(1:length(proposal), function (x) proposal[[x]]$R_sd)

    R_cor.cur=lapply(1:length(current_vals[[1]][[1]]), function (x) current_vals[[1]][[1]][[x]]$R_cor)
    R_sd.cur=lapply(1:length(current_vals[[1]][[1]]), function (x) current_vals[[1]][[1]][[x]]$R_sd)

    hr=sum(unlist(lapply(1:length(current_vals[[1]][[1]]), function (x) proposal[[x]]$hr)))
    jj=sum(unlist(lapply(1:length(current_vals[[1]][[1]]), function (x) proposal[[x]]$jj)))
    cur.liks=current_vals[[1]][[2]]$cur.liks
    #H_fixed=F

    X <- lapply(3:length(pa_data[[1]]), function(j)  pa_data[[1]][[j]])

    #dat.cur.bt=lapply(1:length(dat.cur), function(x) t(apply(dat.cur[[x]], 1, backTransform1)))

    #A.cur.bt=lapply(1:length(A.cur), function(x) backTransform1(A.cur[[x]]))
    #uopdate likelihood efficiently depending on prop

    #print(exists("prop"))
    #print(exists("pa_data"))

    if (prop<3){
      #GLM calc this cou calculating the prior ld be made more efficient for height by
      missSP<-sort(unique(unlist(lapply(1:length(proposal), function(x) proposal[[x]]$missing))))

      missSPdat<-lapply(1:length(missSP), function(x) pa_data[[missSP[x]]])

      #print(pa_data[proposal[[x]]$missing])
      prop_traits_only<-lapply(1:length(A.cur), function(y) dat.prop.bt[[y]][missSP,])
      #curr_traits_missing<-lapply(1:length(A.cur), function(x) dat.cur.bt[[x]][missing,])
      #print(pa_data[proposal[[x]]$missing])
      missing.glm.lik <- GLM_Likelihood_MV_fn(res=prop_traits_only, full_pa_data=missSPdat)
      prop.glm.lik <- cur.liks$glm.lik
      cur.liks$glm.lik[missSP]
      prop.glm.lik[missSP] <- missing.glm.lik
      prop.glm.sum.lik <- sum(unlist(prop.glm.lik))
      GLM <- prop.glm.sum.lik-cur.liks$glm.sum.lik


      if(prop==0){

        #if height was chosen as the prior then you gotta calculate how the prior changes



        prop.prior.height<-lapply(1:length(proposal), function(pred){

          lapply(missSP, function(sp){
            if(is.null(prior[[pred]]$heights.prior)==F){


              prior[[pred]]$heights.prior[[sp]](as.numeric(dat.prop[[pred]][sp,3]))
            }else{
              0
            }}
          )  }
        )





        prior.heights.lik<-cur.liks$prior.lik$prior.heights.lik



        for(pred in 1:length(proposal)){

          names( prior.heights.lik[[pred]])<-tree$tip.label

          prior.heights.lik[[pred]][missSP]<-prop.prior.height[[pred]]


        }


        Prior.height<- sum(unlist(prop.prior.height))-sum(unlist(lapply(1:length(proposal), function(pred) cur.liks$prior.lik$prior.heights.lik[[pred]][missSP])))



      } else{
        #if height wasnt chosen the prior wont change
        prior.heights.lik<-cur.liks$prior.lik$prior.heights.lik
        Prior.height <- 0

      }


    } else {

      prop.glm.lik <- cur.liks$glm.lik

      prop.glm.sum.lik <- cur.liks$glm.sum.lik

      GLM <- 0

      prior.heights.lik<-cur.liks$prior.lik$prior.heights.lik
      Prior.height <- 0

    }


    prop.bm.lik=sum(unlist(lapply(1:length(dat.cur), function(x) lnL_ratematrix(dat.prop[[x]], tree, A.prop[[x]], proposal[[x]]$R))))

    BM <- prop.bm.lik - cur.liks$bm.lik

    prior.root.lik=sum(unlist(lapply(1:length(A.prop.bt), function (x) prior[[x]]$mean.prior(A.prop.bt[[x]]))))
    prior.corr.lik=sum(unlist(lapply(1:length(R_cor.prop), function (x) prior[[x]]$corr.prior(R_cor.prop[[x]]))))
    prior.sd.lik=sum(unlist(lapply(1:length(R_sd.prop), function (x) prior[[x]]$sd.prior(R_sd.prop[[x]]))))

    Prior.mean <-prior.root.lik -  cur.liks$prior.lik$prior.root.lik
    Prior.corr <-prior.corr.lik -  cur.liks$prior.lik$prior.corr.lik
    Prior.sd   <-prior.sd.lik   -  cur.liks$prior.lik$prior.sd.lik

    Prior <- Prior.mean + Prior.corr + Prior.sd + Prior.height

    if (prior_only==T){
      a.ratio=exp(Prior+BM+hr+jj)

    } else if (glm_only==T){
      a.ratio=exp(GLM+hr+jj)

    } else{
      a.ratio=exp(Prior+GLM+BM+hr+jj)
    }

    return(list(a.ratio=a.ratio,
                prop.liks=list(glm.lik=prop.glm.lik,
                               glm.sum.lik=prop.glm.sum.lik,
                               bm.lik=prop.bm.lik,
                               prior.lik=list(prior.root.lik=prior.root.lik, prior.corr.lik=prior.corr.lik, prior.sd.lik=prior.sd.lik, prior.heights.lik=prior.heights.lik)),
                cur.liks=cur.liks,
                hr=hr,
                jj=jj))


  }


  #####Estimator###############################################################################################################################################################################


  metro_haste_full_MV<- function(R_corr_start,
                                 R_sd_start,
                                 A_start,
                                 Prior,
                                 tree,
                                 tibble_data,
                                 pa_data,
                                 iterations,
                                 burnin,
                                 move_prob=c(2/5,2/5,1/15,1/15,1/15),
                                 n,
                                 print.i.freq=100,
                                 print.ac.freq=10,
                                 printing=TRUE,
                                 trim=T,
                                 trim_freq=1,
                                 H_fixed=T,
                                 tuning,
                                 k,
                                 all_props=F,
                                 center_fixed=F,
                                 write_file=T,
                                 IDlen=5,
                                 dir,
                                 outname,
                                 prior_only=F,
                                 glm_only=F,
                                 plot=F,
                                 plot_freq,
                                 plot_file,
                                 True_pars=NULL


  ){
    #set params if line by line#####
    {

      #i=1
      #tree_full<-tree
      #R_corr_start = startPars_scaled$R$R_cor
      #R_sd_start   = startPars_scaled$R$R_sd
      #A_start      = startPars_scaled$A$A_bt
      #Prior=Prior_scale
      #tree=tree
      #tibble_data = startPars_scaled$sim_dat$sim_td_bt
      #pa_data=sets_full$training
      #
      #



      ##beginning real function EVERYTHING BABOVE THIS IS FOR LINE BY LINE

      dat_ft=lapply(1:length(tibble_data), function(pred) t(apply(tibble_data[[pred]], 1, forwardTransform1)))

      #assign tuning
      #
      #
      d<-tuning$niche_prop
      w_mu<-tuning$w_mu
      w_sd<-tuning$w_sd
      v_cor<-tuning$v_cor

      for(pred in 1:length(tibble_data)){

        rownames(dat_ft[[pred]])<-tree$tip.label
      }

      tibble_data_ft<-lapply(1:length(tibble_data), function(x) make.treedata(tree, dat_ft[[x]])$dat)

      A_start_ft<-lapply( 1:length(tibble_data), function(i) forwardTransform1(A_start[[i]])[1:length(A_start[[i]])])

      R<-lapply( 1:length(R_sd_start), function (i) rebuild.cov( r=R_corr_start[[i]], v=(R_sd_start[[i]])^2))



      starting.jacobian <- lapply(1:length(R_sd_start), function(i) sum( sapply(1:length(R_sd_start[i]), function(x) log( ((R_sd_start[[i]][x])^2)) ) ) * log( (length(R_sd_start[[i]])-1)/2 ))
      current_vals=list()
      current_vals[[1]]=list()
      current_vals[[1]][[1]]=list()
      current_vals[[1]][[2]]=list()

      for (i in 1:length(tibble_data)){
        current_vals[[1]][[1]][[i]]=list()
        current_vals[[1]][[1]][[i]]$R=R[[i]]
        current_vals[[1]][[1]][[i]]$R_cor=R_corr_start[[i]]
        current_vals[[1]][[1]][[i]]$R_sd=R_sd_start[[i]]
        current_vals[[1]][[1]][[i]]$A=A_start_ft[[i]]
        current_vals[[1]][[1]][[i]]$dat=tibble_data_ft[[i]]
        current_vals[[1]][[1]][[i]]$curr.jac=starting.jacobian[[i]]
      }

      #likelihood info that appllies across all predictors
      current_vals[[1]][[2]]$move=NULL
      current_vals[[1]][[2]]$accept=NULL
      current_vals[[1]][[2]]$a.ratio=NULL
      current_vals[[1]][[2]]$lik_sum=NULL
      current_vals[[1]][[2]]$cur.liks=NULL

      #keeop track of acceptance ratio proposal through the chain
      current_vals[[1]][[2]]$prop_accept_ratios=NULL

      #makePrior function from ratematrix package produces list of priors, currently most uninformative priors used
      prior <- Prior


      #generate tree.data object for current values (prof add)
      dt=lapply(1:length(current_vals[[1]][[1]]), function(x) data.frame(species=tree$tip.label, as.data.frame(current_vals[[1]][[1]][[x]]$dat)))
      td=lapply(1:length(dt), function(x) make.treedata(tree, dt[[x]]))

      #generate likelihood for current values (prof add)
      #generate likelihood for current values (prof add)

      dat.cur=lapply(1:length(current_vals[[1]][[1]]), function (x) current_vals[[1]][[1]][[x]]$dat)
      dat.cur.bt=lapply(1:length(dat.cur), function(x) t(apply(dat.cur[[x]], 1, backTransform1)))

      A.cur=lapply(1:length(current_vals[[1]][[1]]), function (x) current_vals[[1]][[1]][[x]]$A)
      A.cur.bt=lapply(1:length(A.cur), function(x) backTransform1(A.cur[[x]])[1:length(A.cur[[i]])])


      cur.glm.lik=GLM_Likelihood_MV_fn(dat.cur.bt, pa_data)
      cur.glm.sum.lik=sum(unlist(cur.glm.lik))

      cur.bm.lik=sum(unlist(lapply(1:length(dat.cur), function(x) lnL_ratematrix(dat.cur[[x]], tree, current_vals[[1]][[1]][[x]]$A, current_vals[[1]][[1]][[x]]$R))))


      prior.root.lik=sum(unlist(lapply(1:length(A.cur), function (x) prior[[x]]$mean.prior(A.cur.bt[[x]]))))
      prior.corr.lik=sum(unlist(lapply(1:length(current_vals[[1]][[1]]), function (x) prior[[x]]$corr.prior(current_vals[[1]][[1]][[x]]$R_cor))))
      prior.sd.lik=sum(unlist(lapply(1:length(current_vals[[1]][[1]]), function (x) prior[[x]]$sd.prior(current_vals[[1]][[1]][[x]]$R_sd))))




      prior.heights.lik=lapply(1:length(current_vals[[1]][[1]]), function(pred){
        lapply(1:length(pa_data), function(sp) {
          if(is.null(Prior[[pred]]$heights.prior)==F){
            Prior[[pred]]$heights.prior[[sp]](dat.cur.bt[[pred]][sp,3])
          } else{
            0
          }
        })
      })





      #sum(unlist(lapply(1:length(current_vals[[1]][[1]]), function (x) prior[[x]]$sd.prior(current_vals[[1]][[1]][[x]]$R_sd))))

      #prior[[x]]$sd.prior(c(0.5563417,  0.2075794))

      #c(0.5563417,  0.1075794)

      current_vals[[1]][[2]]$cur.liks=list(glm.lik=cur.glm.lik,
                                           glm.sum.lik=cur.glm.sum.lik,
                                           bm.lik=cur.bm.lik,
                                           prior.lik=list(prior.root.lik=prior.root.lik,
                                                          prior.corr.lik=prior.corr.lik,
                                                          prior.sd.lik=prior.sd.lik,
                                                          prior.heights.lik=prior.heights.lik))

      #keep track of each type of moves accepted, rejected, and rejected with N/A acceptance ratio

      acceptances=0
      A.center.bm.moves=0
      A.center.slide.moves=0
      A.center.mult.moves=0
      A.width.bm.moves=0
      A.width.slide.moves=0
      A.width.mult.moves=0
      A.theta.moves=0
      A.R.corr.moves=0
      A.R.sd.moves=0
      A.R.sd.moves=0

      NA_moves=0
      NA.center.bm.moves=0
      NA.center.slide.moves=0
      NA.center.mult.moves=0
      NA.width.bm.moves=0
      NA.width.slide.moves=0
      NA.width.mult.moves=0
      NA.theta.moves=0
      NA.R.corr.moves=0
      NA.R.sd.moves=0

      Rej.center.bm.moves=0
      Rej.center.slide.moves=0
      Rej.center.mult.moves=0
      Rej.width.bm.moves=0
      Rej.width.slide.moves=0
      Rej.width.mult.moves=0
      Rej.theta.moves=0
      Rej.R.corr.moves=0
      Rej.R.sd.moves=0

      #create objects for chain to be stored in
      chain=list()
      trimmed_chain=list()


      #td_real<-td
      td=lapply(1:length(dt), function(x) make.treedata(tree, dt[[x]]))

      #us proposal list for dbugging to trace throuhg MCMC with the chain
      #proposal_list<-list()
      #set.seed(1)

    }

    ###dev line by line
    #
    #current_vals[[1]]<-results[[3]]$chain[[chain_end]]
    if(all_props==T){
      props<-list()
    }


    ID <- paste( sample(x=1:9, size=IDlen, replace=TRUE), collapse="")

    if(write_file==T){

      Write_to_File(dir=dir,outname=outname,ID=ID,current_vals = current_vals)
    }

    #####startMCMC####
    for (i in 1:iterations){

      k=k


      #   for(w in 1:1){

      for (x in (1:length(dt))){
        td[[x]]$dat<-current_vals[[1]][[1]][[x]]$dat
      }



      prop<-sample(c(0,1,2,3,4,5), size=1, prob = move_prob)


      if(prop<3){

        node =  sample(1:(length(tree$tip.label)+tree$Nnode  ) ,1)
        #extract clade to make shift
        if(node <= length(tree$tip.label)){
          clade <- list(tip.label=tree$tip.label[node])
        } else {
          clade<-extract.clade(tree, node = node )
        }

        clade

        missing<-(1:nrow(td[[x]]$dat))[tree$tip.label%in%clade$tip.label]




      }





      if (prop==0){
        ##center moves
        niche_move=sample(c(1,2), size=1, prob = c(1/2,1/2))

        if(niche_move==1){
          #proposal= lapply(1:length(td), function(x) Slide_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center, niche_move = "center"))
          proposal= lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$slide[missing,3], niche_move = "height", clade))

          #td[[x]]$dat
          #proposal[[x]]$td$dat
          #proposal[[x]]$missing

          move="Height_Slide"

        }else if(niche_move==2){
          #proposal = lapply(1:length(td), function(x) Multiplier_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center, niche_move = "center"))
          proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$multi[missing,3],  niche_move = "height", clade))

          move="Height_Multiplier"
          #lol<-unlist(proposal, recursive=F)
        }


        #there is only one height, need to fix across tips


        #there is only one height, need to fix across tips

        for(z in 2:length(proposal)){

          proposal[[z]]$td$dat[,3]<-proposal[[1]]$td$dat[,3]

          proposal[[z]]$hr=0

        }



      }else if(prop==1){
        ##center moves
        niche_move=sample(c(1,2,3), size=1, prob = c(0,1/2,1/2))

        if(niche_move==1){
          #proposal = lapply(1:length(td), function(x) propMVBM(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, V=NULL))
          move="Center_BM"

        } else if(niche_move==2){
          #proposal= lapply(1:length(td), function(x) Slide_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$center, niche_move = "center"))
          proposal= lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$slide[missing,1] , niche_move = "center", clade))

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
          proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$multi[missing,1], niche_move = "center", clade))
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
          proposal= lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$slide[missing,2], niche_move="width", clade))

          # td[[x]]$dat== proposal[[x]]$td$dat
          #proposal[[x]]$td$dat
          move="Width_Slide"

        }else if(niche_move==3){
          #proposal = lapply(1:length(td), function(x) Multiplier_Proposal_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, H_fixed = H_fixed, d[[x]]))
          proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$multi[missing,2] , niche_move="width", clade))

          move="Width_Multiplier"
          #lol<-unlist(proposal, recursive=F)
        }


      } else if(prop==3){
        ##theta move
        proposal <- lapply(1:length(td), function(x) slideWindow_theta(w_mu[[x]]$slide, tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,  current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac))
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
        proposal <- lapply(1:length(td), function(x) SD_prop(tree, td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac, w_sd=w_sd[[x]]$slide))
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


      probab=acceptance_ratio_MV_fun(proposal,
                                     current_vals,
                                     pa_data,
                                     tree,
                                     prior,
                                     prop,
                                     k,
                                     prior_only=prior_only,
                                     glm_only=glm_only)


      if(is.na(probab$a.ratio)){
        cat(red("NA acceptance ratio ", i+1, " ", current_vals[[1]][[2]]$move, "\n"))
        NA_moves=NA_moves+1

        if (move=="Center_Multiplier"){
          NA.center.mult.moves=NA.center.mult.moves+1
        }else if(move=="Center_Slide"){
          NA.center.slide.moves=NA.center.slide.moves+1
        }else if(move=="Center_BM"){
          NA.center.bm.moves=NA.center.bm.moves+1

        }else if (move=="Width_Multiplier"){
          NA.width.mult.moves=NA.width.mult.moves+1
        }else if(move=="Width_Slide"){
          NA.width.slide.moves=NA.width.slide.moves+1
        }else if(move=="Width_BM"){
          NA.width.bm.moves=NA.width.bm.moves+1

        }else if(move=="theta"){
          NA.theta.moves=NA.theta.moves+1

        }else if(move=="R_corr"){
          NA.R.corr.moves=NA.R.corr.moves+1

        }else if(move=="R_sd"){
          NA.R.sd.moves=NA.R.sd.moves+1
        }

        probab$a.ratio<- -66666666666666
        ##if acceptance ration is above the runif pull it is accepted anc current values are updated
      }
      if (runif(1) < probab$a.ratio){
        for (y in 1:length(current_vals[[1]][[1]])){
          current_vals[[1]][[1]][[y]]$R = proposal[[y]]$R
          current_vals[[1]][[1]][[y]]$R_cor= proposal[[y]]$R_cor
          current_vals[[1]][[1]][[y]]$R_sd= proposal[[y]]$R_sd
          current_vals[[1]][[1]][[y]]$A = proposal[[y]]$A
          current_vals[[1]][[1]][[y]]$dat = proposal[[y]]$td$dat
        }
        current_vals[[1]][[2]]$move = move
        current_vals[[1]][[2]]$accept=TRUE
        current_vals[[1]][[2]]$a.ratio=probab$a.ratio
        current_vals[[1]][[2]]$cur.liks=probab$prop.liks
        current_vals[[1]][[2]]$lik_sum=probab
        current_vals[[1]][[2]]$curr.jac=proposal$prop_jac
        current_vals[[1]][[2]]$prop_accept_ratios= list(center.bm=(A.center.bm.moves/(A.center.bm.moves+Rej.center.bm.moves)),
                                                        center.slide=(A.center.slide.moves/(A.center.slide.moves+Rej.center.slide.moves)),
                                                        center.mult=(A.center.mult.moves/(A.center.mult.moves+Rej.center.mult.moves)),

                                                        width.bm=(A.width.bm.moves/(A.width.bm.moves+Rej.width.bm.moves)),
                                                        width.slide=(A.width.slide.moves/(A.width.slide.moves+Rej.width.slide.moves)),
                                                        width.mult=(A.width.mult.moves/(A.width.mult.moves+Rej.width.mult.moves)),


                                                        theta.slide=(A.theta.moves/(A.theta.moves+Rej.theta.moves)),
                                                        R.corr.slide=(A.R.corr.moves/(A.R.corr.moves+Rej.R.corr.moves)),
                                                        R.sd.slide=(A.R.sd.moves/(A.R.sd.moves+Rej.R.sd.moves)))

        ##keep track fo accepted moves
        acceptances=acceptances+1

        if (move=="Center_Multiplier"){
          A.center.mult.moves=A.center.mult.moves+1
        }else if(move=="Center_Slide"){
          A.center.slide.moves=A.center.slide.moves+1
        }else if(move=="Center_BM"){
          A.center.bm.moves=A.center.bm.moves+1

        }else if (move=="Width_Multiplier"){
          A.width.mult.moves=A.width.mult.moves+1
        }else if(move=="Width_Slide"){
          A.width.slide.moves=A.width.slide.moves+1
        }else if(move=="Width_BM"){
          A.width.bm.moves=A.width.bm.moves+1

        }else if(move=="theta"){
          A.theta.moves=A.theta.moves+1

        }else if(move=="R_corr"){
          A.R.corr.moves=A.R.corr.moves+1

        }else if(move=="R_sd"){
          A.R.sd.moves=A.R.sd.moves+1
        }


        if (acceptances%%print.ac.freq==0 & printing==TRUE){
          #print(proposal[[1]]td$dat)
          cat("iterations:", i+1, "\n")
          cat("acceptances:", acceptances, "\n")
          cat("acceptance ratio:", probab$a.ratio, "\n")
          cat("move utilized:", move, "\n")
        }
        #keep track of rejected moves
      }else {
        #update a.ratio for this step
        current_vals[[1]][[2]]$accept=FALSE
        current_vals[[1]][[2]]$a.ratio=probab$a.ratio
        current_vals[[1]][[2]]$lik_sum=probab
        current_vals[[1]][[2]]$move = move
        current_vals[[1]][[2]]$prop_accept_ratios= list(center.bm=(A.center.bm.moves/(A.center.bm.moves+Rej.center.bm.moves)),
                                                        center.slide=(A.center.slide.moves/(A.center.slide.moves+Rej.center.slide.moves)),
                                                        center.mult=(A.center.mult.moves/(A.center.mult.moves+Rej.center.mult.moves)),

                                                        width.bm=(A.width.bm.moves/(A.width.bm.moves+Rej.width.bm.moves)),
                                                        width.slide=(A.width.slide.moves/(A.width.slide.moves+Rej.width.slide.moves)),
                                                        width.mult=(A.width.mult.moves/(A.width.mult.moves+Rej.width.mult.moves)),


                                                        theta.slide=(A.theta.moves/(A.theta.moves+Rej.theta.moves)),
                                                        R.corr.slide=(A.R.corr.moves/(A.R.corr.moves+Rej.R.corr.moves)),
                                                        R.sd.slide=(A.R.sd.moves/(A.R.sd.moves+Rej.R.sd.moves)))


        if (move=="Center_Multiplier"){
          Rej.center.mult.moves=Rej.center.mult.moves+1
        }else if(move=="Center_Slide"){
          Rej.center.slide.moves=Rej.center.slide.moves+1
        }else if(move=="Center_BM"){
          Rej.center.bm.moves=Rej.center.bm.moves+1

        }else if (move=="Width_Multiplier"){
          Rej.width.mult.moves=Rej.width.mult.moves+1
        }else if(move=="Width_Slide"){
          Rej.width.slide.moves=Rej.width.slide.moves+1
        }else if(move=="Width_BM"){
          Rej.width.bm.moves=Rej.width.bm.moves+1

        }else if(move=="theta"){
          Rej.theta.moves=Rej.theta.moves+1

        }else if(move=="R_corr"){
          Rej.R.corr.moves=Rej.R.corr.moves+1

        }else if(move=="R_sd"){
          Rej.R.sd.moves=Rej.R.sd.moves+1
        }
        if (i%%print.i.freq==0 & printing==TRUE){
          cat("iterations:", i+1, "\n")
          cat("acceptances:", acceptances,"\n")
          cat("acceptance ratio (not accepted):", probab$a.ratio, "\n")
          cat("move utilized:", move, "\n")
        }
      }
      #add current values to chain
      if(i > burnin){
        if(write_file==T){
          if(trim==TRUE){
            if(i%%trim_freq==0){

              Write_to_File(dir=dir, outname=outname, ID=ID, current_vals = current_vals, append=T)

            }
          }else{
            Write_to_File(dir=dir, outname=outname, ID=ID, current_vals = current_vals, append=T)
          }

        }else if(trim==TRUE){
          if(i%%trim_freq==0){
            trimmed_chain=append(trimmed_chain, current_vals[1])

            if(plot==T & i%%plot_freq==0){

              print(length(trimmed_chain))

              pdf(file= paste(plot_file))

              if (is.null(True_pars)==T){
                plotTraces(trimmed_chain=trimmed_chain, plot.true = F)

              }else{
                plotTraces(trimmed_chain=trimmed_chain, plot.true = T, True_pars=True_pars)

              }

              dev.off()
            }



          }
        }else{
          chain=append(chain, current_vals[1])
        }

      }

      if(all_props==T){

        props[[i]]<-list(proposal, move)
      }

      #print(i)

    }


    #length(trimmed_chain)
    ##tail(lapply(1:length(chain), function(x) (chain[[x]][[1]][[1]]$R==chain[[x]][[1]][[2]]$R)))
    ##lapply(1:length(chain), function(x) c(lapply(1:2, function(y) c(chain[[x]][[1]][[y]]$R_sd, chain[[x]][[1]][[y]]$A)), chain[[x]][[2]]$a.ratio, chain[[x]][[2]]$lik_sum$prop.liks$prior.lik  ) )
    #lapply(1:2, function(y)  trimmed_chain[[4500]][[1]][[y]]$R_sd)
    #
    #TruePars_scale[[1]]$A
    #
    #plot(unlist(unique(lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[2]]))))
    #
    #plot(unlist(unique(lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$A))))
    #
    ##lapply(1:2, function(y) chain[[x]][[1]][[y]]$A)
    ##
    ##lapply(1:length(chain), function(x) c(chain[[x]][[2]]$move, chain[[x]][[2]]$accept))
    ##lapply(1:length(chain), function(x) chain[[x]][[2]]$accept)
    #dev.off()
    #chain[[1]][[1]][[1]]$A

    if(all_props==T){
      return(
        list(accept_ratios = list(center.bm=(A.center.bm.moves/(A.center.bm.moves+Rej.center.bm.moves)),
                                  center.slide=(A.center.slide.moves/(A.center.slide.moves+Rej.center.slide.moves)),
                                  center.mult=(A.center.mult.moves/(A.center.mult.moves+Rej.center.mult.moves)),

                                  width.bm=(A.width.bm.moves/(A.width.bm.moves+Rej.width.bm.moves)),
                                  width.slide=(A.width.slide.moves/(A.width.slide.moves+Rej.width.slide.moves)),
                                  width.mult=(A.width.mult.moves/(A.width.mult.moves+Rej.width.mult.moves)),


                                  theta.slide=(A.theta.moves/(A.theta.moves+Rej.theta.moves)),
                                  R.corr.slide=(A.R.corr.moves/(A.R.corr.moves+Rej.R.corr.moves)),
                                  R.sd.slide=(A.R.sd.moves/(A.R.sd.moves+Rej.R.sd.moves))),

             accepted.moves=list(A.center.bm.moves   = A.center.bm.moves,
                                 A.center.slide.moves= A.center.slide.moves,
                                 A.center.mult.moves = A.center.mult.moves,

                                 A.width.bm.moves    = A.width.bm.moves,
                                 A.width.slide.moves = A.width.slide.moves,
                                 A.width.mult.moves  = A.width.mult.moves,

                                 A.theta.moves       = A.theta.moves,
                                 A.R.corr.moves      = A.R.corr.moves,
                                 A.R.sd.moves        = A.R.sd.moves),

             rejected.moves=list(Rej.center.bm.moves   = Rej.center.bm.moves,
                                 Rej.center.slide.moves= Rej.center.slide.moves,
                                 Rej.center.mult.moves = Rej.center.mult.moves,

                                 Rej.width.bm.moves    = Rej.width.bm.moves,
                                 Rej.width.slide.moves = Rej.width.slide.moves,
                                 Rej.width.mult.moves  = Rej.width.mult.moves,

                                 Rej.theta.moves       = Rej.theta.moves,
                                 Rej.R.corr.moves      = Rej.R.corr.moves,
                                 Rej.R.sd.moves        = Rej.R.sd.moves),

             chain=trimmed_chain,

             NAs=NA_moves,

             NA.moves= list(NA.center.bm.moves  = NA.center.bm.moves,
                            NA.center.slide.moves= NA.center.slide.moves,
                            NA.center.mult.moves = NA.center.mult.moves,

                            NA.width.bm.moves    = NA.width.bm.moves,
                            NA.width.slide.moves = NA.width.slide.moves,
                            NA.width.mult.moves  = NA.width.mult.moves,

                            NA.theta.moves       = NA.theta.moves,
                            NA.R.corr.moves      = NA.R.corr.moves,
                            NA.R.sd.moves        = NA.R.sd.moves),

             all_props=props,

             "full")
      )

    }



    if(trim==TRUE){
      return(list(acceptances=acceptances,

                  accept_ratios = list(center.bm=(A.center.bm.moves/(A.center.bm.moves+Rej.center.bm.moves)),
                                       center.slide=(A.center.slide.moves/(A.center.slide.moves+Rej.center.slide.moves)),
                                       center.mult=(A.center.mult.moves/(A.center.mult.moves+Rej.center.mult.moves)),

                                       width.bm=(A.width.bm.moves/(A.width.bm.moves+Rej.width.bm.moves)),
                                       width.slide=(A.width.slide.moves/(A.width.slide.moves+Rej.width.slide.moves)),
                                       width.mult=(A.width.mult.moves/(A.width.mult.moves+Rej.width.mult.moves)),


                                       theta.slide=(A.theta.moves/(A.theta.moves+Rej.theta.moves)),
                                       R.corr.slide=(A.R.corr.moves/(A.R.corr.moves+Rej.R.corr.moves)),
                                       R.sd.slide=(A.R.sd.moves/(A.R.sd.moves+Rej.R.sd.moves))),

                  accepted.moves=list(A.center.bm.moves   = A.center.bm.moves,
                                      A.center.slide.moves= A.center.slide.moves,
                                      A.center.mult.moves = A.center.mult.moves,

                                      A.width.bm.moves    = A.width.bm.moves,
                                      A.width.slide.moves = A.width.slide.moves,
                                      A.width.mult.moves  = A.width.mult.moves,

                                      A.theta.moves       = A.theta.moves,
                                      A.R.corr.moves      = A.R.corr.moves,
                                      A.R.sd.moves        = A.R.sd.moves),

                  rejected.moves=list(Rej.center.bm.moves   = Rej.center.bm.moves,
                                      Rej.center.slide.moves= Rej.center.slide.moves,
                                      Rej.center.mult.moves = Rej.center.mult.moves,

                                      Rej.width.bm.moves    = Rej.width.bm.moves,
                                      Rej.width.slide.moves = Rej.width.slide.moves,
                                      Rej.width.mult.moves  = Rej.width.mult.moves,

                                      Rej.theta.moves       = Rej.theta.moves,
                                      Rej.R.corr.moves      = Rej.R.corr.moves,
                                      Rej.R.sd.moves        = Rej.R.sd.moves),

                  chain=trimmed_chain,

                  NAs=NA_moves,

                  NA.moves= list(NA.center.bm.moves  = NA.center.bm.moves,
                                 NA.center.slide.moves= NA.center.slide.moves,
                                 NA.center.mult.moves = NA.center.mult.moves,

                                 NA.width.bm.moves    = NA.width.bm.moves,
                                 NA.width.slide.moves = NA.width.slide.moves,
                                 NA.width.mult.moves  = NA.width.mult.moves,

                                 NA.theta.moves       = NA.theta.moves,
                                 NA.R.corr.moves      = NA.R.corr.moves,
                                 NA.R.sd.moves        = NA.R.sd.moves),

                  "full"))
    }
    else{


      return(list(acceptances=acceptances,

                  accept_ratios = list(center.bm=(A.center.bm.moves/(A.center.bm.moves+Rej.center.bm.moves)),
                                       center.slide=(A.center.slide.moves/(A.center.slide.moves+Rej.center.slide.moves)),
                                       center.mult=(A.center.mult.moves/(A.center.mult.moves+Rej.center.mult.moves)),

                                       width.bm=(A.width.bm.moves/(A.width.bm.moves+Rej.width.bm.moves)),
                                       width.slide=(A.width.slide.moves/(A.width.slide.moves+Rej.width.slide.moves)),
                                       width.mult=(A.width.mult.moves/(A.width.mult.moves+Rej.width.mult.moves)),

                                       theta.slide=(A.theta.moves/(A.theta.moves+Rej.theta.moves)),
                                       R.corr.slide=(A.R.corr.moves/(A.R.corr.moves+Rej.R.corr.moves)),
                                       R.sd.slide=(A.R.sd.moves/(A.R.sd.moves+Rej.R.sd.moves))),

                  accepted.moves=list(A.center.bm.moves   = A.center.bm.moves,
                                      A.center.slide.moves= A.center.slide.moves,
                                      A.center.mult.moves = A.center.mult.moves,

                                      A.width.bm.moves    = A.width.bm.moves,
                                      A.width.slide.moves = A.width.slide.moves,
                                      A.width.mult.moves  = A.width.mult.moves,

                                      A.theta.moves       = A.theta.moves,
                                      A.R.corr.moves      = A.R.corr.moves,
                                      A.R.sd.moves        = A.R.sd.moves),

                  rejected.moves=list(Rej.center.bm.moves   = Rej.center.bm.moves,
                                      Rej.center.slide.moves= Rej.center.slide.moves,
                                      Rej.center.mult.moves = Rej.center.mult.moves,
                                      Rej.width.bm.moves    = Rej.width.bm.moves,
                                      Rej.width.slide.moves = Rej.width.slide.moves,
                                      Rej.width.mult.moves  = Rej.width.mult.moves,
                                      Rej.theta.moves       = Rej.theta.moves,
                                      Rej.R.corr.moves      = Rej.R.corr.moves,
                                      Rej.R.sd.moves        = Rej.R.sd.moves),

                  chain=chain,

                  NAs=NA_moves,

                  NA.moves=list(NA.center.bm.moves   = NA.center.bm.moves,
                                NA.center.slide.moves= NA.center.slide.moves,
                                NA.center.mult.moves = NA.center.mult.moves,
                                NA.width.bm.moves    = NA.width.bm.moves,
                                NA.width.slide.moves = NA.width.slide.moves,
                                NA.width.mult.moves  = NA.width.mult.moves,
                                NA.theta.moves       = NA.theta.moves,
                                NA.R.corr.moves      = NA.R.corr.moves,
                                NA.R.sd.moves        = NA.R.sd.moves),

                  "full"))

    }
  }





  ##MCMC output################################################################################################################################################################################################################################################################################################################################################

  MCMC_df_fn<- function(MCMC_results, FUN="median"){
    all.dat<-   lapply(1:length(MCMC_results$chain[[1]][[1]]), function(pred) lapply(1:length(MCMC_results$chain), function(x)   t(apply(MCMC_results$chain[[x]][[1]][[pred]]$dat,1,backTransform1))))
    listVec <- lapply(1:length(MCMC_results$chain[[1]][[1]]), function(pred) lapply(all.dat[[pred]], c, recursive=TRUE))
    m <- lapply(1:length(MCMC_results$chain[[1]][[1]]), function(pred) do.call(cbind, listVec[[pred]]) )
    out <- lapply(1:length(MCMC_results$chain[[1]][[1]]), function(pred) as.data.frame(matrix(apply(m[[pred]], 1, FUN), nrow=nrow(MCMC_results$chain[[1]][[1]][[1]]$dat), ncol=ncol(MCMC_results$chain[[1]][[1]][[1]]$dat))))
    return(out)

  }




  chain2mcmcobj_full_H_fixed_MV<- function(pred, chain, tree, object, zeroed=FALSE, true_data=NULL, true_A=NULL, true_R_sd=NULL, true_R_cor=NULL){
    #object specifies if converting a trait from the response curve, theta, or the R matrix
    #zeroing is an option only if using simulated data
    i=pred
    if (object=="trait"){
      if (zeroed==TRUE){
        tmp=lapply(1:length(chain[[5]]), function(x) t(chain[[5]][[x]][[1]][[i]]$dat-true_data[[i]]))
      } else{
        tmp <- lapply(1:length(chain[[5]]), function(x) t(apply(chain[[5]][[x]][[1]][[i]]$dat,1,backTransform1)))
      }
      tmp2 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x))))
      mymcmc_1 <- mcmc(tmp2[,1:length(tree$tip.label)])
      mymcmc_2 <-mcmc(tmp2[,(length(tree$tip.label)+1):(2*length(tree$tip.label))])
      mymcmc_3 <-mcmc(tmp2[,(2*length(tree$tip.label)+1):(3*length(tree$tip.label))])
      return(list(mymcmc.center=mymcmc_1, mymcmc.width=mymcmc_2, mymcmc.height=mymcmc_3))
    } else if( object=="theta"){
      if (zeroed==TRUE){
        tmp <- lapply(1:length(chain[[5]]), function(x) backTransform1(chain[[5]][[x]][[1]][[i]]$A)-true_A[[i]])
      }else{
        tmp <-  lapply( 1:length(chain[[5]]),  function(x) backTransform1(chain[[5]][[x]][[1]][[i]]$A))
      }
      A_1 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[1]]))))
      A_2 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[2]]))))
      A_3 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[3]]))))
      mymcmc_A_1<-mcmc(A_1)
      mymcmc_A_2<-mcmc(A_2)
      mymcmc_A_3<-mcmc(A_3)
      return(list(mymcmc_A.center=mymcmc_A_1, mymcmc_A.width=mymcmc_A_2, mymcmc_A.height=mymcmc_A_3))

    } else if( object=="R_sd"){
      if (zeroed==TRUE){
        tmp <-  lapply( 1:length(chain[[5]]),  function(x) chain[[5]][[x]][[1]][[i]]$R_sd-true_R_sd[[i]])
      }else{
        tmp <-  lapply( 1:length(chain[[5]]),  function(x) chain[[5]][[x]][[1]][[i]]$R_sd)
      }
      sd_1 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[1]]))))
      sd_2 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[2]]))))
      #sd_3 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[3]]))))
      mymcmc_sd_1<-mcmc(sd_1)
      mymcmc_sd_2<-mcmc(sd_2)
      #mymcmc_sd_3<-mcmc(sd_3)
      return(list(mymcmc_sd.center=mymcmc_sd_1, mymcmc_sd.width=mymcmc_sd_2))
    } else if( object=="R_corr"){
      if (zeroed==TRUE){
        tmp <- lapply(chain[[5]], function(x) x$R_cor-true_R[[i]])
      } else{
        tmp <- lapply( 1:length(chain[[5]]),  function(x) chain[[5]][[x]][[1]][[i]]$R_cor)
      }
      tmp_1_2<-lapply(tmp, function(x) x[[1,2]])
      #tmp_1_3<-lapply(tmp, function(x) x[[1,3]])
      #tmp_2_3<-lapply(tmp, function(x) x[[2,3]])
      corr_1_2 <-  do.call(rbind, lapply(tmp_1_2, function(x) unlist(data.frame(x))))
      #corr_1_3 <-  do.call(rbind, lapply(tmp_1_3, function(x) unlist(data.frame(x))))
      #corr_2_3 <-  do.call(rbind, lapply(tmp_2_3, function(x) unlist(data.frame(x))))
      mymcmc_1_2<-mcmc(corr_1_2)
      #mymcmc_1_3<-mcmc(corr_1_3)
      #mymcmc_2_3<-mcmc(corr_2_3)

      return(list(mymcmc_c.w=mymcmc_1_2))
    }

    return()
  }





  BePhyNE_out2coda_mcmc<- function(mcmc_output, tree){
    #iterates through all predictors and saves chains as coda mcmc objects that are easy to plot



    mcmc_object<-list()
    mcmc_object$tip_curves<-list()
    mcmc_object$A<-list()
    mcmc_object$R_sd<-list()
    mcmc_object$R_cor<-list()



    nenvir_pred<-length(mcmc_output$chain[[1]][[1]])


    for(pred in 1: nenvir_pred){

      #response curve mcmc chains

      tmp <- lapply(1:length(mcmc_output$chain), function(x) t(apply(mcmc_output$chain[[x]][[1]][[pred]]$dat,1,backTransform1)))

      tmp2 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x))))
      mymcmc_1 <- mcmc(tmp2[,1:length(tree$tip.label)])
      mymcmc_2 <-mcmc(tmp2[,(length(tree$tip.label)+1):(2*length(tree$tip.label))])
      mymcmc_3 <-try(mcmc(tmp2[,(2*length(tree$tip.label)+1):(3*length(tree$tip.label))]),silent = T)
      mcmc_object$tip_curves[[pred]]<-list(optimum=mymcmc_1, breadth=mymcmc_2, tolerance=mymcmc_3)


      # mvBM root means

      tmp <-  lapply( 1:length(mcmc_output$chain),  function(x) backTransform1((mcmc_output$chain[[x]][[1]][[pred]]$A)))
      A_1 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[1]]))))
      A_2 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[2]]))))
      A_3 <-  try(do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[3]])))),silent = T)
      mymcmc_A_1<-mcmc(A_1)
      mymcmc_A_2<-mcmc(A_2)
      mymcmc_A_3<-try(mcmc(A_3),silent = T)
      mcmc_object$A[[pred]]<-list(optimum=mymcmc_A_1, breadth=mymcmc_A_2, tolerance=mymcmc_A_3)




      # mvBM rates
      tmp <-  lapply( 1:length(mcmc_output$chain),  function(x) mcmc_output$chain[[x]][[1]][[pred]]$R_sd)
      sd_1 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[1]]))))
      sd_2 <-  do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[2]]))))
      sd_3 <-  try(do.call(rbind, lapply(tmp, function(x) unlist(data.frame(x[[3]])))),silent = T)
      mymcmc_sd_1<-mcmc(sd_1)
      mymcmc_sd_2<-mcmc(sd_2)
      mymcmc_sd_3<-try(mcmc(sd_3),silent = T)
      mcmc_object$R_sd[[pred]]<-list(optimum=mymcmc_sd_1, breadth=mymcmc_sd_2, tolerance=mymcmc_sd_3 )



      # mvBM Correlation matrix


      tmp <-  lapply( 1:length(mcmc_output$chain),  function(x)(mcmc_output$chain[[x]][[1]][[pred]]$R_cor))
      tmp_1_2<-lapply(tmp, function(x) x[[1,2]])
      tmp_1_3<-try(lapply(tmp, function(x) x[[1,3]]),silent = T)
      tmp_2_3<-try(lapply(tmp, function(x) x[[2,3]]),silent = T)
      cor_1_2 <-  do.call(rbind, lapply(tmp_1_2, function(x) unlist(data.frame(x))))
      cor_1_3 <-  try(do.call(rbind, lapply(tmp_1_3, function(x) unlist(data.frame(x)))),silent = T)
      cor_2_3 <-  try(do.call(rbind, lapply(tmp_2_3, function(x) unlist(data.frame(x)))),silent = T)
      mymcmc_1_2<-mcmc(cor_1_2)
      mymcmc_1_3<-try(mcmc(cor_1_3),silent = T)
      mymcmc_2_3<-try(mcmc(cor_2_3),silent = T)

      mcmc_object$R_cor[[pred]]<-list(opt.bre=mymcmc_1_2, opt.tol=  mymcmc_1_3, br_tol=mymcmc_2_3)

    }


    names(mcmc_object$tip_curves)<- unlist(lapply(1:nenvir_pred, function(pred) paste("pred", pred, sep="_")))
    names(mcmc_object$A)<- lapply(1:nenvir_pred, function(pred) paste("pred", pred, sep="_"))
    names(mcmc_object$R_sd)<- lapply(1:nenvir_pred, function(pred) paste("pred", pred, sep="_"))
    names(mcmc_object$R_cor)<- lapply(1:nenvir_pred, function(pred) paste("pred", pred, sep="_"))

    return(mcmc_object)
  }






  HPD_list<-function(mcmc_object, tree, prob=0.95){

    #creates list of HPDs each list entry is named after the specific parameter the HPD describes

    hpd_list<-lapply(unlist(unlist(mcmc_lt, recursive = F), recursive = F), function(c) try(HPDinterval(c, args=list(prob=prob)),silent = T))


    for(df in grep("tip_curves",x = names(hpd_list))){

      rownames(hpd_list[[df]])<-tree$tip.label

    }

    hpd_list_final<-hpd_list[!lapply(hpd_list,class)=="try-error"]

    return(hpd_list_final)

  }




  ##Predict#######################################################################################################################################################################################################################

  ##predict and cross validation functions

  '%!in%' <- function(x,y)!('%in%'(x,y))



  separate.data<-function(data_final, ratio=.5){

    training_set=list()
    testing_set=list()

    for( sp in 1:length(data_final)){


      if(is.na(data_final[[sp]]$y[[1]])==T){

        training_set[[sp]]<-data_final[[sp]]
        testing_set[[sp]] <-data_final[[sp]]



      }else{

        df<-data.frame(data_final[[sp]])

        pres<-which(df$y==1)




        training_pres<-sample(pres, length(pres)*ratio)
        if(length(training_pres)==0){
          training_pres<-sample(pres, 1)
        }

        testing_pres<-pres[pres%!in%training_pres]


        abs<-which(df$y==0)
        training_abs<-sample(abs, length(abs)*ratio)
        testing_abs<-abs[abs%!in%training_abs]


        training_set[[sp]]<-as.list(df[sort(c(training_abs, training_pres)),])
        training_set[[sp]]$species<-as.character(training_set[[sp]]$species[[1]])

        testing_set[[sp]]<-as.list(df[sort(c(testing_abs, testing_pres)),])
        testing_set[[sp]]$species<-as.character(testing_set[[sp]]$species[[1]])



      }

    }

    return(list(training=training_set, predicting=testing_set))

  }




  predict.ENE<-function(traits, pa_data){

    #must use pa_data like object with list of list format with data for each sp

    X= lapply(1:nrow(traits[[1]]), function(sp)  list(pa_data[[sp]]$X1, pa_data[[sp]]$X2))

    betas<-lapply(1:length(traits), function(pred) traits2coefs(traits[[pred]]))

    ###glm stuff
    yy_prelogit_sep  = lapply( 1:nrow(traits[[1]]), function(sp)  lapply(1:length(traits), function(pred) (betas[[pred]][sp,1] + betas[[pred]][sp,2]*X[[sp]][[pred]]+ betas[[pred]][sp,3]*(X[[sp]][[pred]]^2)  ) ) )

    yy_prelogit= lapply(1:nrow(traits[[1]]), function(sp) rowSums(matrix(unlist(yy_prelogit_sep[[sp]]), ncol=length(yy_prelogit_sep[[sp]]), byrow=F)))

    presProb <-lapply(1:nrow(traits[[1]]), function(sp) data.frame(y=pa_data[[sp]]$y, fitted.values=1/(1+exp(-1*yy_prelogit[[sp]] ) )) ) #convert to probability using logit link
    ###

    return(presProb)

  }





  assess.predict<-function(presProb, plot=F){

    if(is.na(presProb$y)==T){

      paste("no data for this one chief")
      return(list(AUC=NA, TSS=NA, ctable=NA))
    }

    eval.input=presProb
    cutoffs=seq(0,1,.02)
    cutoffs=as.matrix(cutoffs)
    classify=apply(cutoffs, 1, function(y) apply(eval.input, 1, function(x) if(x[1]== 1 && x[2] >= y) {'TP'} else if (x[1]== 0 && x[2] < y){'TN'} else if (x[1]== 0 && x[2] >= y){'FP'} else {'FN'}))
    ctable=apply(classify, 2, function(x) c(length(which(x == "TP")),length(which(x == "TN")),length(which(x == "FP")),length(which(x == "FN"))))
    #View(ctable)
    ctable=t(ctable)
    confusions=cbind(cutoffs,ctable)
    colnames(confusions)=c("Cutoff","Correct Event (TP)","Correct Non-Event (TN)","Incorrect Event (FP)", "Incorrect Non-Event (FN)")

    ###Formulas- defns: Correct Event (TP); Correct Non-Event (TN); Incorrect Event (FP); Incorrect Non-Event (FN)
    ##FCC=(TP+TN)/n
    ##Sensitivity=TP/(TP+FN)
    ##Specificity=TN/(TN+FP)
    ##False Positive Rate =FP/(TP+FP)
    ##False Negative Rate =FN/(TN+FN)

    FCC=(confusions[,2]+confusions[,3])/nrow(eval.input)*100;Sensitivity=confusions[,2]/(confusions[,2]+confusions[,5])*100;Specificity=confusions[,3]/(confusions[,3]+confusions[,4])*100;False_POS_rate=confusions[,4]/(confusions[,2]+confusions[,4])*100;False_NEG_rate=confusions[,5]/(confusions[,5]+confusions[,3])*100
    roc_x=100-Specificity
    ctable_complete=cbind(confusions,FCC,Sensitivity,Specificity,False_POS_rate,False_NEG_rate,roc_x)



    #install.packages("flux")
    library(flux)
    AUC=auc(roc_x, Sensitivity, thresh = NULL, dens = 100)/100

    TSS=Sensitivity+Specificity-100

    if (plot==T){

      plot(roc_x,Sensitivity,type="l")
      abline(0,1)


      plot(cutoffs, FCC, type='p', col="white", xlab="Probability Cutoff", ylab="Rate(%)", xlim=c(0,1), ylim=c(0,100))
      lines(cutoffs,FCC, lty=1, col="blue",lwd=2)
      lines(cutoffs,Sensitivity, lty=2, col="red",lwd=2)
      lines(cutoffs,Specificity, lty=3, col="black",lwd=2)
      lines(cutoffs,TSS, lty=4, col="green",lwd=2)
      legend(x="topright", c("Fraction Correctly Classified (FCC)","Sensitivity","Specificity",
                             "True Skill Statistic (TSS)"),lty=1:3,col=c("blue","red","black", "green"), cex=1)
      abline(v=0.38, col="gold", lwd=4)
    }

    return(list(AUC=AUC, TSS=TSS, ctable=ctable_complete))

  }


  ####Plotting##############################################################################################################################################################################################################################################################################################################################################

  plotTraitsContour<-function(X1_unique, X2_unique,  traits){


    X1=as.vector(matrix(X1_unique, nrow = length(X1_unique) ,ncol = length(X1_unique), byrow = T))
    X2=as.vector(matrix(X2_unique, nrow = length(X2_unique) ,ncol = length(X2_unique), byrow = F))


    yyp<-list()
    p_mat<-list()

    #backtransformed trait values
    res <-traits
    betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients

    plot(0, 0, type="n", xlim=c(min(X1), max(X1)), ylim=c(min(X2), max(X2)), xlab="Pred_1", ylab="Pred_2")


    for (i in 1:nrow(res[[1]])){

      yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*X1 + betas[[1]][i,3]*X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*X2 + betas[[2]][i,3]*X2^2)

      yyp[[i]] <- 1/(1+exp(-1*yy)) #convert to probability

      p_mat[[i]]<-matrix(yyp[[i]], nrow=length(X1_unique), ncol=length(X2_unique), byrow=T)


      pred=p_mat[[i]]
      #pred=exp(PL[[i]])/(1+exp(PL[[i]]))
      max(pred)
      min(pred)

      #data.interp <- interp(X1,X2,pred)
      #image(data.interp) #or
      #contour(data.interp, col=i, add=T)

      #plot(X, cex=pred)
      contour(unique(X1), unique(X2), pred, col=i, add=T)
      points(res[[1]][[i,1]], res[[2]][[i,1]], pch=21, cex=2, bg= i)
      arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, col="red", lty=2)
      arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, col="red", lty=2)

    }
  }



  PlotContourLik<-function(pa_data, surface_traits, contours=NULL, surface_traitscont_col="black", contour_cols=c('red', 'green'), fixed_presence=T , rep=NULL){


    #i=1
    ##
    ##X1_unique=full_X1
    ##X2_unique=full_X2
    ###
    #surface_traits =list(TruePars_scale[[i]]$sim_dat$sim_dat_bt[[1]],
    #                   TruePars_scale[[i]]$sim_dat$sim_dat_bt[[2]])
    #surface_traits=true
    #contours=list(start = start[[i]]
    #              ,final =list(t(apply(results[[i]]$chain[[chain_end]][[1]][[1]]$dat, 1, backTransform1)),
    #                          t(apply(results[[i]]$chain[[chain_end]][[1]][[2]]$dat, 1, backTransform1)))
    #)
    #
    ##contour_cols<-c("red", "green")
    ##
    ##
    #pa_data=pa_data_scaled[[i]]$species_data
    #
    #
    #pa_data[[i]]$X1
    #pa_data[[i]]$X2

    #yyp<-list()
    #p_mat<-list()

    #likelihood gets plotted for prop
    res <- surface_traits
    res_contours <- contours

    betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
    betas_contours <- lapply(1:length(res_contours), function(i) lapply(1:length(res_contours[[i]]), function(x) traits2coefs( res_contours[[i]][[x]]) ) )# Convert to beta coefficients


    if(fixed_presence==T){

      if(is.null(contours)==T){

        for (i in 1:nrow(res[[1]])){

          yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1 + betas[[1]][i,3]*pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2 + betas[[2]][i,3]*pa_data[[i]]$X2^2)

          yyp <- 1/(1+exp(-1*yy)) #convert to probability

          pl= dbinom(1,1,prob= yyp, log=T)

          yyp_mat<-matrix(yyp, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)


          ##commented out unless for line by line checking pred matrix if something seems fucked
          #dimnames(yyp_mat) <- list(unique(pa_data[[i]]$X1),unique(pa_data[[i]]$X2) )
          #
          #View(yyp_mat)

          p_mat<-matrix(pl, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

          yyp_pred=yyp_mat
          lik_pred=p_mat
          #pred=exp(PL[[i]])/(1+exp(PL[[i]]))
          #max(pred)
          #min(pred)

          #data.interp <- interp(X1,X2,pred)
          #image(data.interp) #or
          #contour(data.interp, col=i, add=T)

          #plot(X, cex=pred)

          #contour(unique(X1), unique(X2), lik_pred,  add=T)


          filled.contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), lik_pred ,
                         plot.title = {
                           Arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                           Arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                           contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred, col=surface_traitscont_col,  add=T);
                           title(main = paste("likelihood surface fixed pres", "rep", rep, "species", i))
                         })

        }
      }else{

        for (i in 1:nrow(res[[1]])){


          yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1 + betas[[1]][i,3]*pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2 + betas[[2]][i,3]*pa_data[[i]]$X2^2)

          yyp <- 1/(1+exp(-1*yy)) #convert to probability

          pl= dbinom(1,1,prob= yyp, log=T)

          yyp_mat<-matrix(yyp, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=T)

          p_mat<-matrix(pl, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=T)

          yy_contours = lapply(1:length(res_contours), function(x)  (betas_contours[[x]][[1]][i,1] + betas_contours[[x]][[1]][i,2]*pa_data[[i]]$X1 + betas_contours[[x]][[1]][i,3]*pa_data[[i]]$X1^2) + (betas_contours[[x]][[2]][i,1] + betas_contours[[x]][[2]][i,2]*pa_data[[i]]$X2 + betas_contours[[x]][[2]][i,3]*pa_data[[i]]$X2^2) )

          yyp_contours<- lapply(1:length(res_contours), function(x) 1/(1+exp(-1*yy_contours[[x]])) )#convert to probability

          yyp_mat_contours<-lapply(1:length(res_contours), function(x) matrix(yyp_contours[[x]], nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=T) )


          yyp_pred=yyp_mat
          lik_pred=p_mat
          #pred=exp(PL[[i]])/(1+exp(PL[[i]]))
          #max(pred)
          #min(pred)

          #data.interp <- interp(X1,X2,pred)
          #image(data.interp) #or
          #contour(data.interp, col=i, add=T)

          #plot(X, cex=pred)

          #contour(unique(X1), unique(X2), lik_pred,  add=T)
          filled.contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), lik_pred ,
                         plot.title = {
                           Arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                           Arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                           contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred, col=surface_traitscont_col,  add=T);

                           lapply(1:length(res_contours), function(x) Arrows(res_contours[[x]][[1]][[i,1]] +res_contours[[x]][[1]][[i,2]], res_contours[[x]][[2]][[i,1]], res_contours[[x]][[1]][[i,1]] - res_contours[[x]][[1]][[i,2]], res_contours[[x]][[2]][[i,1]], code = 3, arr.length =.1  ,col=contour_cols[[x]], lty=1, lwd=2) );
                           lapply(1:length(res_contours), function(x) Arrows(res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]] - res_contours[[x]][[2]][[i,2]], res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]] + res_contours[[x]][[2]][[i,2]], code = 3, arr.length =.1  ,col=contour_cols[[x]], lty=1, lwd=2) );
                           lapply(1:length(res_contours), function(x) contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_mat_contours[[x]], col=contour_cols[[x]],  add=T) );

                           title(main = paste("likelihood surface fixed pres", "rep", rep, "species", i))
                         })


          #contour(unique(X1)-.33, unique(X2), yyp_pred, col=i,  add=T)
          #plot(0, 0, type="n", xlim=c(min(X1), max(X1)), ylim=c(min(X2), max(X2)), xlab="Pred_1", ylab="Pred_2")

          #filled.contour(unique(X1), unique(X2), pred,  add=T)
          #points(res[[1]][[i,1]]-.33, res[[2]][[i,1]], pch=21, cex=2, bg= i)
          #Arrows(res[[1]][[i,1]] +res[[1]][[i,2]]-.33, res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]]-.33, res[[2]][[i,1]], code = 3, arr.length =.1  ,col="black", lty=1, lwd=2)
          #Arrows(res[[1]][[i,1]]-.33, res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]]-.33, res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col="black", lty=1, lwd=2)
        }

      }
    }else{
      if(is.null(contours)==T){
        for (i in 1:nrow(res[[1]])){

          yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1 + betas[[1]][i,3]*pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2 + betas[[2]][i,3]*pa_data[[i]]$X2^2)

          yyp <- 1/(1+exp(-1*yy)) #convert to probability

          yyp_mat<-matrix(yyp, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

          pl= dbinom(pa_data[[i]]$y, 1, prob = yyp, log=T)

          p_mat<-matrix(pl, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

          yyp_pred=yyp_mat

          lik_pred=p_mat

          filled.contour( unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), lik_pred ,
                          plot.title = {
                            Arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                            Arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                            contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred, col=1,  add=T);

                            title(main = paste("likelihood surface", "rep", rep, "species", i))

                          }
          )
        }

      }else{

        for (i in 1:nrow(res[[1]])){


          yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1 + betas[[1]][i,3]*pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2 + betas[[2]][i,3]*pa_data[[i]]$X2^2)

          yyp <- 1/(1+exp(-1*yy)) #convert to probability


          yyp_mat<-matrix(yyp, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)


          #length(unique(pa_data[[i]]$X1))
          #length(unique(pa_data[[i]]$X2))

          yy_contours = lapply(1:length(res_contours), function(x)  (betas_contours[[x]][[1]][i,1] + betas_contours[[x]][[1]][i,2]*pa_data[[i]]$X1 + betas_contours[[x]][[1]][i,3]*pa_data[[i]]$X1^2) + (betas_contours[[x]][[2]][i,1] + betas_contours[[x]][[2]][i,2]*pa_data[[i]]$X2 + betas_contours[[x]][[2]][i,3]*pa_data[[i]]$X2^2) )

          yyp_contours<- lapply(1:length(res_contours), function(x) 1/(1+exp(-1*yy_contours[[x]])) )#convert to probability

          yyp_mat_contours<-lapply(1:length(res_contours), function(x) matrix(yyp_contours[[x]], nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F) )


          pl= dbinom(pa_data[[i]]$y, 1, prob = yyp, log=T)

          p_mat<-matrix(pl, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)


          yyp_pred=yyp_mat

          lik_pred=p_mat
          #pred=exp(PL[[i]])/(1+exp(PL[[i]]))

          max(yyp_pred)
          min(yyp_pred)

          max(lik_pred)
          min(lik_pred)

          #data.interp <- interp(X1,X2,pred)
          #image(data.interp) #or
          #contour(data.interp, col=i, add=T)

          #plot(X, cex=pred)

          #contour(unique(X1), unique(X2), lik_pred,  add=T)
          filled.contour( unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), lik_pred ,
                          plot.title = {
                            Arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                            Arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col=surface_traitscont_col, lty=1, lwd=2);
                            contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred, col=1,  add=T);

                            lapply(1:length(res_contours), function(x) Arrows(res_contours[[x]][[1]][[i,1]] +res_contours[[x]][[1]][[i,2]], res_contours[[x]][[2]][[i,1]], res_contours[[x]][[1]][[i,1]] - res_contours[[x]][[1]][[i,2]], res_contours[[x]][[2]][[i,1]], code = 3, arr.length =.1  ,col=contour_cols[[x]], lty=1, lwd=2) );
                            lapply(1:length(res_contours), function(x) Arrows(res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]] - res_contours[[x]][[2]][[i,2]], res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]] + res_contours[[x]][[2]][[i,2]], code = 3, arr.length =.1  ,col=contour_cols[[x]], lty=1, lwd=2) );
                            lapply(1:length(res_contours), function(x) contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_mat_contours[[x]], col=contour_cols[[x]],  add=T) );

                            title(main = paste("likelihood surface", "rep", rep, "species", i))

                          }
          )

          #contour(unique(X1)-.33, unique(X2), yyp_pred, col=i,  add=T)
          #plot(0, 0, type="n", xlim=c(min(X1), max(X1)), ylim=c(min(X2), max(X2)), xlab="Pred_1", ylab="Pred_2")

          #filled.contour(unique(X1), unique(X2), pred,  add=T)
          #points(res[[1]][[i,1]]-.33, res[[2]][[i,1]], pch=21, cex=2, bg= i)
          #Arrows(res[[1]][[i,1]] +res[[1]][[i,2]]-.33, res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]]-.33, res[[2]][[i,1]], code = 3, arr.length =.1  ,col="black", lty=1, lwd=2)
          #Arrows(res[[1]][[i,1]]-.33, res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]]-.33, res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col="black", lty=1, lwd=2)
        }
      }
    }
  }

  PlotContourPred<-function(pa_data, surface_traits, contours=NULL,  cont_col="black", contour_cols=c('red', 'green'), rep=NULL){


    #i=1
    ##
    ##X1_unique=full_X1
    ##X2_unique=full_X2
    ###
    #surface_traits =list(TruePars_scale[[i]]$sim_dat$sim_dat_bt[[1]],
    #                   TruePars_scale[[i]]$sim_dat$sim_dat_bt[[2]])
    ###
    #contours=list(start = start[[i]]
    #final =list(t(apply(results[[i]]$chain[[chain_end]][[1]][[1]]$dat, 1, backTransform1)),
    #           t(apply(results[[i]]$chain[[chain_end]][[1]][[2]]$dat, 1, backTransform1)))
    #)
    ##
    #contour_cols<-c("red", "green")
    #
    #cont_col="black"
    ##
    #  ##
    #i=1
    #
    #  pa_data=pa_data_scaled[[i]]$species_data
    #
    #  surface_traits =true[[i]]
    #  #pa_data[[i]]$X1
    #pa_data[[i]]$X2

    #yyp<-list()
    #p_mat<-list()





    #likelihood gets plotted for prop
    res <- surface_traits
    res_contours <- contours

    betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
    betas_contours <- lapply(1:length(res_contours), function(i) lapply(1:length(res_contours[[i]]), function(x) traits2coefs( res_contours[[i]][[x]]) ) )# Convert to beta coefficients


    if(is.null(contours)==T){

      for (i in 1:nrow(res[[1]])){




        yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1 + betas[[1]][i,3]*pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2 + betas[[2]][i,3]*pa_data[[i]]$X2^2)

        yyp <- 1/(1+exp(-1*yy)) #convert to probability

        yyp_mat<-matrix(yyp, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

        pl= dbinom(pa_data[[i]]$y, 1, prob = yyp, log=T)

        p_mat<-matrix(pl, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

        yyp_pred=yyp_mat

        lik_pred=p_mat

        filled.contour( unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred ,
                        plot.title = {points(res[[1]][[i,1]], res[[2]][[i,1]], pch=21, cex=2, bg= i);
                          Arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, arr.length =.1  ,col=cont_col, lty=1, lwd=2);
                          Arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col=cont_col, lty=1, lwd=2);
                          contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred, col=1,  add=T);

                          title(main = paste("likelihood surface", "rep", rep, "species", i))

                        }
        )
      }

    }else{

      for (i in 1:nrow(res[[1]])){


        yy  =  (betas[[1]][i,1] + betas[[1]][i,2]*pa_data[[i]]$X1 + betas[[1]][i,3]*pa_data[[i]]$X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*pa_data[[i]]$X2 + betas[[2]][i,3]*pa_data[[i]]$X2^2)

        yyp <- 1/(1+exp(-1*yy)) #convert to probability


        yyp_mat<-matrix(yyp, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

        #length(unique(pa_data[[i]]$X1))
        #length(unique(pa_data[[i]]$X2))

        yy_contours = lapply(1:length(res_contours), function(x)  (betas_contours[[x]][[1]][i,1] + betas_contours[[x]][[1]][i,2]*pa_data[[i]]$X1 + betas_contours[[x]][[1]][i,3]*pa_data[[i]]$X1^2) + (betas_contours[[x]][[2]][i,1] + betas_contours[[x]][[2]][i,2]*pa_data[[i]]$X2 + betas_contours[[x]][[2]][i,3]*pa_data[[i]]$X2^2) )

        yyp_contours<- lapply(1:length(res_contours), function(x) 1/(1+exp(-1*yy_contours[[x]])) )#convert to probability

        yyp_mat_contours<-lapply(1:length(res_contours), function(x) matrix(yyp_contours[[x]], nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F) )

        pl= dbinom(pa_data[[i]]$y, 1, prob = yyp, log=T)

        p_mat<-matrix(pl, nrow=length(unique(pa_data[[i]]$X1)), ncol=length(unique(pa_data[[i]]$X2)), byrow=F)

        yyp_pred=yyp_mat

        lik_pred=p_mat
        #pred=exp(PL[[i]])/(1+exp(PL[[i]]))

        max(yyp_pred)
        min(yyp_pred)

        max(lik_pred)
        min(lik_pred)

        #data.interp <- interp(X1,X2,pred)
        #image(data.interp) #or
        #contour(data.interp, col=i, add=T)

        #plot(X, cex=pred)

        #contour(unique(X1), unique(X2), lik_pred,  add=T)
        filled.contour( unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred ,
                        plot.title = {points(res[[1]][[i,1]], res[[2]][[i,1]], pch=21, cex=2, bg= i);
                          Arrows(res[[1]][[i,1]] +res[[1]][[i,2]], res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]], res[[2]][[i,1]], code = 3, arr.length =.1  ,col=cont_col, lty=1, lwd=2);
                          Arrows(res[[1]][[i,1]], res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]], res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col=cont_col, lty=1, lwd=2);
                          contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_pred, col=1,  add=T);

                          lapply(1:length(res_contours), function(x) points(res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]], pch=21, cex=2, bg= i) );
                          lapply(1:length(res_contours), function(x) Arrows(res_contours[[x]][[1]][[i,1]] +res_contours[[x]][[1]][[i,2]], res_contours[[x]][[2]][[i,1]], res_contours[[x]][[1]][[i,1]] - res_contours[[x]][[1]][[i,2]], res_contours[[x]][[2]][[i,1]], code = 3, arr.length =.1  ,col=contour_cols[[x]], lty=1, lwd=2) );
                          lapply(1:length(res_contours), function(x) Arrows(res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]] - res_contours[[x]][[2]][[i,2]], res_contours[[x]][[1]][[i,1]], res_contours[[x]][[2]][[i,1]] + res_contours[[x]][[2]][[i,2]], code = 3, arr.length =.1  ,col=contour_cols[[x]], lty=1, lwd=2) );
                          lapply(1:length(res_contours), function(x) contour(unique(pa_data[[i]]$X1), unique(pa_data[[i]]$X2), yyp_mat_contours[[x]], col=contour_cols[[x]],  add=T) );

                          title(main = paste("likelihood surface", "rep", rep, "species", i))

                        }
        )

        #contour(unique(X1)-.33, unique(X2), yyp_pred, col=i,  add=T)
        #plot(0, 0, type="n", xlim=c(min(X1), max(X1)), ylim=c(min(X2), max(X2)), xlab="Pred_1", ylab="Pred_2")

        #filled.contour(unique(X1), unique(X2), pred,  add=T)
        #points(res[[1]][[i,1]]-.33, res[[2]][[i,1]], pch=21, cex=2, bg= i)
        #Arrows(res[[1]][[i,1]] +res[[1]][[i,2]]-.33, res[[2]][[i,1]], res[[1]][[i,1]] - res[[1]][[i,2]]-.33, res[[2]][[i,1]], code = 3, arr.length =.1  ,col="black", lty=1, lwd=2)
        #Arrows(res[[1]][[i,1]]-.33, res[[2]][[i,1]] - res[[2]][[i,2]], res[[1]][[i,1]]-.33, res[[2]][[i,1]] + res[[2]][[i,2]], code = 3, arr.length =.1  ,col="black", lty=1, lwd=2)
      }
    }

  }

  CompPlotTraitsRidge<-function(tree, X1_unique, X2_unique,  traits, traits2, traits3){
    TL <- max(branching.times(tree)) # Get total tree height

    maxX <- 2.5*TL # Set max X limit for plot
    xplot <- c(1.5*TL, 2.5*TL) #Range of where the curves will go
    #res <- lapply(1:length(results$trimmed_chain), function(x) backTransform1(results$trimmed_chain[[x]]$dat)) #Backtransform the posterior distribution
    #betas <- traits2coefs(res[[gen]]) # Convert to beta coefficients

    X1=X1_unique
    X2=X2_unique



    X<-list(X1, X2)

    yyp<-list()
    p_mat<-list()

    #backtransformed trait values
    res <-traits
    res2<-traits2
    res3<-traits3
    betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
    betas2<-  lapply(1:length(res), function(x) traits2coefs(res2[[x]]))
    betas3<-  lapply(1:length(res), function(x) traits2coefs(res3[[x]]))

    # trait 1 is red
    # trait 2 is green
    # trait 3 is black

    yy<-list()
    yy[[1]]  = lapply(1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*X1 + betas[[1]][i,3]*X1^2) + (betas[[2]][i,1] + betas[[2]][i,2]*res[[2]][[i,1]]+ betas[[2]][i,3]*res[[2]][[i,1]]^2))
    yy[[2]] = lapply(1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*res[[1]][[i,1]] + betas[[1]][i,3]*res[[1]][[i,1]]^2) + (betas[[2]][i,1] + betas[[2]][i,2]*X2 + betas[[2]][i,3]*X2^2))

    yy2<-list()
    yy2[[1]]  = lapply(1:nrow(res2[[1]]), function(i) (betas2[[1]][i,1] + betas2[[1]][i,2]*X1 + betas2[[1]][i,3]*X1^2) + (betas2[[2]][i,1] + betas2[[2]][i,2]*res2[[2]][[i,1]]+ betas2[[2]][i,3]*res2[[2]][[i,1]]^2))
    yy2[[2]] = lapply(1:nrow(res2[[1]]), function(i) (betas2[[1]][i,1] + betas2[[1]][i,2]*res2[[1]][[i,1]] + betas2[[1]][i,3]*res2[[1]][[i,1]]^2) + (betas2[[2]][i,1] + betas2[[2]][i,2]*X2 + betas2[[2]][i,3]*X2^2))

    yy3<-list()
    yy3[[1]]  = lapply(1:nrow(res3[[1]]), function(i) (betas3[[1]][i,1] + betas3[[1]][i,2]*X1 + betas3[[1]][i,3]*X1^2) + (betas3[[2]][i,1] + betas3[[2]][i,2]*res3[[2]][[i,1]]+ betas3[[2]][i,3]*res3[[2]][[i,1]]^2))
    yy3[[2]] = lapply(1:nrow(res3[[1]]), function(i) (betas3[[1]][i,1] + betas3[[1]][i,2]*res3[[1]][[i,1]] + betas3[[1]][i,3]*res3[[1]][[i,1]]^2) + (betas3[[2]][i,1] + betas3[[2]][i,2]*X2 + betas3[[2]][i,3]*X2^2))


    for(x in 1:length(res)){


      plot(tree, x.lim=c(0,maxX), y.lim=c(0,length(tree$tip.label)+3)) #Plot the tree, with extra X limits for curves and y limits for text
      text(0, length(tree$tip.label)+2, pos=4, cex=1, labels=paste0("Pred = ",  (1:length(res))[[x]] )) #Give the log likelihood
      #text(0, length(tree$tip.label)+2, pos=4, cex=2, labels=paste0("LnLik = ",  round(results$trimmed_chain[[gen]]$cur.liks$glm.sum.lik+results$trimmed_chain[[gen]]$cur.liks$bm.lik),2)) #Give the log likelihood


      for (i in 1:nrow(res[[1]])){


        yyp <- 1/(1+exp(-1*yy[[x]][[i]])) #convert to probability

        xxp <- seq(xplot[1], xplot[2], length.out=length(yyp)) #Shift the true X's (in Celsius) over to the plot XX's, which need to be right of the tree
        lines(xxp, yyp+i, col=3) #Plot the lines
        temp <- X[[x]] # Get the data temperature values
        #pa <- data[[i]]$y #get the presence/absence values from the data
        temp.p <- ( temp - min(X[[x]]) ) /diff(range(X[[x]])) * diff(xplot) + min(xplot) #Rescale the temperature values to appear in the right spot to the right of the tree
        #points(temp.p, 0.5+rep(i, length(pa)), pch="|", col=pa+1, cex=0.2) # Plot the data and color by presence/absence

        yyp2 <- 1/(1+exp(-1*yy2[[x]][[i]])) #convert to probability

        xxp2 <- seq(xplot[1], xplot[2], length.out=(yyp2)) #Shift the true X's (in Celsius) over to the plot XX's, which need to be right of the tree
        lines(xxp2, yyp2+i, col=2) #Plot the lines
        temp2 <- X[[x]] # Get the data temperature values
        #pa <- data[[i]]$y #get the presence/absence values from the data
        #temp.p2 <- ( temp - min(X[[x]]) ) /diff(range(X[[x]])) * diff(xplot) + min(xplot) #Rescale the temperature values to appear in the right spot to the right of the tree
        #points(temp.p, 0.5+rep(i, length(pa)), pch="|", col=pa+1, cex=0.2) # Plot the data and color by presence/absence



        yyp3 <- 1/(1+exp(-1*yy3[[x]][[i]])) #convert to probability

        xxp3 <- seq(xplot[1], xplot[2], length.out=(yyp3)) #Shift the true X's (in Celsius) over to the plot XX's, which need to be right of the tree
        lines(xxp3, yyp3+i, cex=1.25) #Plot the lines
        temp3 <- X[[x]] # Get the data temperature values
        #pa <- data[[i]]$y #get the presence/absence values from the data
        #temp.p3 <- ( temp3 - min(X[[x]]) ) /diff(range(X[[x]])) * diff(xplot) + min(xplot) #Rescale the temperature values to appear in the right spot to the right of the tree
        #points(temp.p, 0.5+rep(i, length(pa)), pch="|", col=pa+1, cex=0.2) # Plot the data and color by presence/absence


      }

    }

  }

  plotContourPA_dat<-function(tree, pa_data, traits, cont_colors, unscale=F, scale_atr=NULL){

    #pa_data<-data_final
    #traits=list(start[[1]],median)
    #tree=phylo
    #

    X1_unique<-seq(-4,4,.01)
    X2_unique<-seq(-4,4,.01)

    if (unscale==T){
      X1_unique<- X1_unique*scale_atr$scale[[1]]+scale_atr$center[[1]]
      X2_unique<- X2_unique*scale_atr$scale[[2]]+scale_atr$center[[2]]
      traits<- lapply(1:length(traits), function(set) unscale(traits[[set]], scale_atr) )

      for (sp in 1:length(pa_data)){
        pa_data[[sp]]$X1<-pa_data[[sp]]$X1*scale_atr$scale[[1]]+scale_atr$center[[1]]
        pa_data[[sp]]$X2<-pa_data[[sp]]$X2*scale_atr$scale[[2]]+scale_atr$center[[2]]
      }
    }



    X1=as.vector(matrix(X1_unique, nrow = length(X1_unique) ,ncol = length(X1_unique), byrow = T))
    X2=as.vector(matrix(X2_unique, nrow = length(X2_unique) ,ncol = length(X2_unique), byrow = F))


    yyp<-list()
    yyp_pred=list()
    p_mat<-list()

    #backtransformed trait values
    for (set in 1:length(traits)){
      res <-traits[[set]]
      betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients



      #likelihood gets plotted for prop

      yy  = lapply(1:length(tree$tip.label), function(sp) (betas[[1]][sp,1] + betas[[1]][sp,2]*X1 + betas[[1]][sp,3]*X1^2) + (betas[[2]][sp,1] + betas[[2]][sp,2]*X2 + betas[[2]][sp,3]*X2^2))

      yyp <-  lapply(1:length(tree$tip.label), function(sp) 1/(1+exp(-1*yy[[sp]])) ) #convert to probability

      yyp_mat<- lapply(1:length(tree$tip.label), function(sp) matrix(yyp[[sp]], nrow=length( X1_unique), ncol=length( X2_unique), byrow=T))

      #pl= lapply(1:length(phylo$tip.label), function(sp) rbinom(size=1, prob = yyp[[sp]]) )

      #p_mat<-lapply(1:length(phylo$tip.label), function(sp) matrix(pl[[sp]], nrow=length(X1_unique), ncol=length(X2_unique), byrow=T) )

      yyp_pred[[set]]=yyp_mat


    }

    for (sp in 1:length(tree$tip.label)){
      plot(pa_data[[sp]]$X1[pa_data[[sp]]$y==0],pa_data[[sp]]$X2[pa_data[[sp]]$y==0],col=2, xlim=c(min(X1_unique)-1,max(X1_unique)+1), ylim=c(min(X2_unique)-1,max(X2_unique)+1), xlab="precip cm", ylab="temp C", main=paste(tree$tip.label[[sp]]))
      points(pa_data[[sp]]$X1[pa_data[[sp]]$y==1],pa_data[[sp]]$X2[pa_data[[sp]]$y==1],col=3)

      lapply(1:length(traits), function(set) contour(X1_unique, X2_unique, yyp_pred[[set]][[sp]], col=cont_colors[[set]],  add=T))

    }


    return()

  }

  plotTraces<-function(trimmed_chain, plot.true=T, True_pars){

    if(plot.true==T){

      par(mfrow = c(3, 4))

      true=mean(True_pars$sim_dat$sim_dat_bt[[1]][,1])


      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,1])))), type="b", cex=.5, main= paste( "center_1 avg true=", "\n", round(true,5) ), xlab="iteration", ylab="center mean" )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,1])))) ), type="b", cex=.5,  main="center_1_avg dens", xlab="center_avg", ylab="dens" )
      abline(v=true, col=3)



      true=True_pars$R$R_sd[[1]][[1]]

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[1]]))) , type="b", cex=.5,  main=paste("Rsd_center_1 true=", "\n",round(true, 5)), xlab="iteration", ylab="center_rsd" )
      abline(h=true, col=3)
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[1]]))) ), type="b", cex=.5,  main="Rsd_center_1 dens", xlab="Rsd_center_1", ylab="dens" )
      abline(v=true,col=3)

      true=True_pars$A$A_bt[[1]][[1]]

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[1]))), type="b", cex=.5,  main=paste("A_center_1 true=", "\n",round(true, 5)), xlab="iteration", ylab="center_A"  )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[1]))) ), type="b", cex=.5,  main="A_center_1 dens", xlab="A_center_1", ylab="dens"  )
      abline(v=true, col=3)

      true=mean(True_pars$sim_dat$sim_dat_bt[[1]][,2])

      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,2])))), type="b", cex=.5, main=paste("width_1 avg true=", "\n",round(true, 5)), xlab="iteration", ylab="width mean" )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,2])))) ), type="b", cex=.5,  main="width_1_avg dens", xlab="width_avg", ylab="dens" )
      abline(v=true, col=3)

      true=True_pars$R$R_sd[[1]][[2]]
      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[2]]))) , type="b", cex=.5,  main=paste("Rsd_width_1 true=", "\n",round(true, 5)), xlab="iteration", ylab="width_rsd" )
      abline(h=true, col=3)
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[2]]))) ), type="b", cex=.5,  main="Rsd_width_1 dens", xlab="Rsd_width_1", ylab="dens" )
      abline(v=true, col=3)


      true=True_pars$A$A_bt[[1]][[2]]
      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[2]))), type="b", cex=.5,  main=paste("A_width_1  true=", "\n",round(true, 5)), xlab="iteration", ylab="width_A"  )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[2]))) ), type="b", cex=.5,  main="A_width_1 dens", xlab="A_width_1", ylab="dens"  )
      abline(v=true, col=3)





      ###pred 2

      true=mean(True_pars$sim_dat$sim_dat_bt[[2]][,1])


      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,1])))), type="b", cex=.5, main= paste( "center_2 avg true=", "\n", round(true,5) ), xlab="iteration", ylab="center mean" )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,1])))) ), type="b", cex=.5,  main="center_2_avg dens", xlab="center_avg", ylab="dens" )
      abline(v=true, col=3)



      true=True_pars$R$R_sd[[2]][[1]]

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[1]]))) , type="b", cex=.5,  main=paste("Rsd_center_2 true=", "\n",round(true, 5)), xlab="iteration", ylab="center_rsd" )
      abline(h=true, col=3)
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[1]]))) ), type="b", cex=.5,  main="Rsd_center_2 dens", xlab="Rsd_center_2", ylab="dens" )
      abline(v=true,col=3)

      true=True_pars$A$A_bt[[2]][[1]]

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[2]]$A)[1]))), type="b", cex=.5,  main=paste("A_center_2 true=", "\n",round(true, 5)), xlab="iteration", ylab="center_A"  )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[2]]$A)[1]))) ), type="b", cex=.5,  main="A_center_2 dens", xlab="A_center_2", ylab="dens"  )
      abline(v=true, col=3)

      true=mean(True_pars$sim_dat$sim_dat_bt[[2]][,2])

      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,2])))), type="b", cex=.5, main=paste("width_2 avg true=", "\n",round(true, 5)), xlab="iteration", ylab="width mean" )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,2])))) ), type="b", cex=.5,  main="width_2_avg dens", xlab="width_avg", ylab="dens" )
      abline(v=true, col=3)

      true=True_pars$R$R_sd[[2]][[2]]
      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[2]]))) , type="b", cex=.5,  main=paste("Rsd_width_2 true=", "\n",round(true, 5)), xlab="iteration", ylab="width_rsd" )
      abline(h=true, col=3)
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[2]]))) ), type="b", cex=.5,  main="Rsd_width_2 dens", xlab="Rsd_width_2", ylab="dens" )
      abline(v=true, col=3)


      true=True_pars$A$A_bt[[2]][[2]]
      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[2]]$A)[2]))), type="b", cex=.5,  main=paste("A_width_2  true=", "\n",round(true, 5)), xlab="iteration", ylab="width_A"  )
      abline(h=true, col=3)
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[2]]$A)[2]))) ), type="b", cex=.5,  main="A_width_2 dens", xlab="A_width_2", ylab="dens"  )
      abline(v=true, col=3)

    } else {

      par(mfrow = c(3, 4))
      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,1])))), type="b", cex=.5, main="center_1 avg", xlab="iteration", ylab="center mean" )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,1])))) ), type="b", cex=.5,  main="center_1_avg dens", xlab="center_avg", ylab="dens" )

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[1]]))) , type="b", cex=.5,  main="Rsd_center_1", xlab="iteration", ylab="center_rsd" )
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[1]]))) ), type="b", cex=.5,  main="Rsd_center_1 dens", xlab="Rsd_center_1", ylab="dens" )

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[1]))), type="b", cex=.5,  main="A_center_1", xlab="iteration", ylab="center_A"  )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[1]))) ), type="b", cex=.5,  main="A_center_1 dens", xlab="A_center_1", ylab="dens"  )

      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,2])))), type="b", cex=.5, main="width_1 avg", xlab="iteration", ylab="width mean" )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[1]]$dat, 1, backTransform1) )[,2])))) ), type="b", cex=.5,  main="width_1_avg dens", xlab="width_avg", ylab="dens" )

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[2]]))) , type="b", cex=.5,  main="Rsd_width_1", xlab="iteration", ylab="width_rsd" )
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[1]]$R_sd[[2]]))) ), type="b", cex=.5,  main="Rsd_width_1 dens", xlab="Rsd_width_1", ylab="dens" )

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[2]))), type="b", cex=.5,  main="A_width_1", xlab="iteration", ylab="width_A"  )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[2]))) ), type="b", cex=.5,  main="A_width_1 dens", xlab="A_width_1", ylab="dens"  )





      ###pred 2

      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,1])))), type="b", cex=.5, main="center_2 avg", xlab="iteration", ylab="center mean" )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,1])))) ), type="b", cex=.5,  main="center_2_avg dens", xlab="center_avg", ylab="dens" )

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[1]]))) , type="b", cex=.5,  main="Rsd_center_2", xlab="iteration", ylab="center_rsd" )
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[1]]))) ), type="b", cex=.5,  main="Rsd_center_2 dens", xlab="Rsd_center_2", ylab="dens" )


      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[2]]$A)[1]))), type="b", cex=.5,  main="A_center_2", xlab="iteration", ylab="center_A"  )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[2]]$A)[1]))) ), type="b", cex=.5,  main="A_center_2 dens", xlab="A_center_2", ylab="dens"  )


      plot( unlist((lapply(1:length(trimmed_chain), function(y) mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,2])))), type="b", cex=.5, main="width_2 avg", xlab="iteration", ylab="width mean" )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)   mean(t(apply(trimmed_chain[[y]][[1]][[2]]$dat, 1, backTransform1) )[,2])))) ), type="b", cex=.5,  main="width_2_avg dens", xlab="width_avg", ylab="dens" )

      plot( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[2]]))) , type="b", cex=.5,  main="Rsd_width_2", xlab="iteration", ylab="width_rsd" )
      plot( density( unlist((lapply(1:length(trimmed_chain), function(y)  trimmed_chain[[y]][[1]][[2]]$R_sd[[2]]))) ), type="b", cex=.5,  main="Rsd_width_2 dens", xlab="Rsd_width_2", ylab="dens" )


      plot( unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[2]))), type="b", cex=.5,  main="A_width_2", xlab="iteration", ylab="width_A"  )
      plot( density(unlist((lapply(1:length(trimmed_chain), function(y)  backTransform1(trimmed_chain[[y]][[1]][[1]]$A)[2]))) ), type="b", cex=.5,  main="A_width_2 dens", xlab="A_width_2", ylab="dens"  )


    }

  }




  PlotBePhyNEchains<-function(mcmc_object, dir=NA, name="BePhyNE", traces=T, densities=T){


    if(is.na(dir)==T){

      dir= paste(getwd(),"/", sep="")
    }

    nenvir_pred<-length(mcmc_object$tip_curves)

    ####traces
    if(traces==T){


      #response curve traits

      {
        pdf(file= paste(dir, name,"curve_traces.pdf",sep="_"))

        par(mfrow = c(3, 2), oma=c(0,0,2,0))


        for(x in 1:length(tree$tip.label)){

          for(i in 1:nenvir_pred){

            traceplot( mcmc_object$tip_curves[[i]]$optimum[,x], ask=F, main=paste("optimum", i ) )

            traceplot( mcmc_object$tip_curves[[i]]$breadth[,x], ask=F, main=paste("breadth", i ) )

            traceplot( mcmc_object$tip_curves[[i]]$tolerance[,x], ask=F, main=paste("tolerance", i ) )


          }
          title( tree$tip.label[[x]] ,outer = T)

        }





        dev.off()
      }

      #A traces

      {

        pdf(file= paste(dir, name,"A_traces.pdf", sep="_"))


        for(i in 1:nenvir_pred){
          par(mfrow = c(3, 2))


          traceplot( mcmc_object$A[[i]]$optimum, ask=F, main=paste("optimum", i ) )

          traceplot( mcmc_object$A[[i]]$breadth, ask=F, main=paste("breadth", i ) )

          #try(traceplot( mcmc_object$A[[i]]$tolerance, ask=F, main=paste("tolerance", i ) ), silent = T)



        }

        dev.off()
      }

      #R_sd rtraces
      {

        pdf(file= paste(dir, name,"R_sd_traces.pdf", sep="_"))


        for(i in 1:nenvir_pred){
          par(mfrow = c(3, 2))


          traceplot( mcmc_object$R_sd[[i]]$optimum, ask=F, main=paste("optimum", i ) )

          traceplot( mcmc_object$R_sd[[i]]$breadth, ask=F, main=paste("breadth", i ) )

          #try(traceplot( mcmc_object$R_sd[[i]]$tolerance, ask=F, main=paste("tolerance", i ) ), silent = T)



        }

        dev.off()
      }

      ##R correlations traces

      {

        pdf(file= paste(dir, name,"R_cor_traces.pdf", sep="_"))


        for(i in 1:nenvir_pred){
          par(mfrow = c(3, 2))


          traceplot( mcmc_object$R_cor[[i]]$opt.bre, ask=F, main=paste("optimum and breadth", i ) )

          # try(traceplot( mcmc_object$R_cor[[i]]$opt.tol, ask=F, main=paste("optimum and tolerance", i ) ), silent = T)

          #try(traceplot( mcmc_object$R_cor[[i]]$br_tol, ask=F, main=paste("breadth and tolerance", i ) ), silent = T)



        }

        dev.off()
      }




      ####end traces
    }



    ######densities

    if(densities==T){

      {
        pdf(file= paste(dir, name,"curve_dens.pdf", sep="_"))

        for(x in 1:length(tree$tip.label)){


          for(i in 1:nenvir_pred){
            par(mfrow = c(2, 2))
            densplot( mcmc_object$tip_curves[[i]]$optimum[,x], ask=F, main=paste("optimum", i ) )

            densplot( mcmc_object$tip_curves[[i]]$breadth[,x], ask=F, main=paste("breadth", i ) )

            densplot( mcmc_object$tip_curves[[i]]$tolerance[,x], ask=F, main=paste("tolerance", i ) )

            title( tree$tip.label[[x]] )

          }

          title( tree$tip.label[[x]] ,outer = T)


        }
        dev.off()
      }

      {

        pdf(file= paste(dir, name,"A_dens.pdf", sep="_"))


        for(i in 1:nenvir_pred){
          par(mfrow = c(3, 2))
          densplot( mcmc_object$A[[i]]$optimum, ask=F, main=paste("optimum", i ) )

          densplot( mcmc_object$A[[i]]$breadth, ask=F, main=paste("breadth", i ) )

          #try(densplot( mcmc_object$A[[i]]$tolerance, ask=F, main=paste("tolerance", i ) ), silent = T)


        }

        dev.off()
      }

      #R_sd rtraces
      {

        pdf(file= paste(dir, name,"R_sd_dens.pdf", sep="_"))


        for(i in 1:nenvir_pred){
          par(mfrow = c(3, 2))


          densplot( mcmc_object$R_sd[[i]]$optimum, ask=F, main=paste("optimum", i ) )

          densplot( mcmc_object$R_sd[[i]]$breadth, ask=F, main=paste("breadth", i ) )

          #try(densplot( mcmc_object$R_sd[[i]]$tolerance, ask=F, main=paste("tolerance", i ) ), silent = T)



        }

        dev.off()
      }

      ##R correlations traces

      {

        pdf(file= paste(dir, name,"R_cor_dens.pdf", sep="_"))


        for(i in 1:nenvir_pred){
          par(mfrow = c(3, 2))


          densplot( mcmc_object$R_cor[[i]]$opt.bre, ask=F, main=paste("optimum and breadth", i ) )

          #try(densplot( mcmc_object$R_cor[[i]]$opt.tol, ask=F, main=paste("optimum and tolerance", i ) ), silent = T)

          #try(densplot( mcmc_object$R_cor[[i]]$br_tol, ask=F, main=paste("breadth and tolerance", i ) ), silent = T)



        }

        dev.off()


      }



    }

  }


  #########file writing ##############
  Write_to_File<- function(dir, outname, ID, current_vals, append=F){

    ## Open files to write:
    files <- list( file(file.path(dir, paste(outname,".",ID,".mcmc",sep="")), open="a"),
                   file(file.path(dir, paste(outname,".",ID,".log",sep="")), open="a")
    )

    if (append==F){
      names<-lapply(1:length(current_vals[[1]][[1]]), function(x) names(unlist(current_vals[[1]][[1]][[x]])) )
      header<-unlist(lapply(1:length(names), function(trait) lapply(1:length(names[[1]]), function(x) paste("pred_" , trait, "_", names[[trait]][[x]], sep=""))))
      cat(header, sep=" ", file=files[[1]], append=TRUE) ## Write the header to the file.
      cat("\n", file=files[[1]], append=TRUE)


      header <- names(unlist(current_vals[[1]][[2]]))
      cat(header, sep=" ", file=files[[2]], append=TRUE) ## Write the header to the file.
      cat("\n", file=files[[2]], append=TRUE)
    }

    cat( c(unname(unlist(current_vals[[1]][[1]]))), sep=" ", file=files[[1]], append=TRUE)
    cat(" ", sep="", file=files[[1]], append=TRUE)
    cat("\n", file=files[[1]], append=TRUE)


    cat(c(unname(unlist(current_vals[[1]][[2]]))), sep=" ", file=files[[2]], append=TRUE)
    cat(" ", sep="", file=files[[2]], append=TRUE)
    cat("\n", file=files[[2]], append=TRUE)

  }


  readNicheMCMC <- function(outname, burn = 0.1, thin = 1, dir=NULL){

    if(is.null(dir)==F){
      setwd(dir)
    }

    ## In this version the posterior is in a single file. iterate through all ID #'s in directory
    mcmc <- read_lines( file=file.path( paste(outname, ".", ID, ".mcmc", sep="")) )

    header <- mcmc[1]
    header <- as.character( strsplit(x=header, split=" ", fixed=TRUE)[[1]] )

    ## Apply thinning and burnin.
    obs.gen <- length( mcmc )-1 ## Compute observed number of samples (Fix for unfinished chains.) Note that header is first line.
    ## Check if we can read the MCMC. In cases when MCMC failed or something.
    if( obs.gen <= thin ) stop("Length of mcmc samples is lower than 'thin'.")
    if( burn <= 0 ){
      post <- seq(2, obs.gen+1, by=thin) ## First line is the header.
    } else{
      post <- seq(round(obs.gen * burn)+1, obs.gen+1, by=thin) ## First line is the header.
    }
    mcmc <- mcmc[post]

    #View(mcmc)

    ## Parse the posterior samples.
    mcmc_new <- mcmc(do.call(rbind, lapply(mcmc, function(x) as.numeric( strsplit(x=x, split=" ", fixed=TRUE)[[1]] )
    )))


    #names<-lapply(1:length(current_vals[[1]][[1]]), function(x) names(unlist(current_vals[[1]][[1]][[x]])) )
    #
    #header<-unlist(lapply(1:length(names), function(trait) lapply(1:length(names[[1]]), function(x) paste("pred_" , trait, "_", names[[trait]][[x]], sep=""))))


    colnames(mcmc_new)<-header

    traits<-lapply(1:length(current_vals[[1]][[1]]), function(trait) mcmc_new[ , grepl( paste("pred_", trait, "_dat.", sep="") , colnames( mcmc_new ) ) ] )

    R_sd<-lapply(1:length(current_vals[[1]][[1]]), function(trait) mcmc_new[ , grepl( paste("pred_", trait, "_R_sd.", sep="") , colnames( mcmc_new ) ) ] )

    R_cor<-lapply(1:length(current_vals[[1]][[1]]), function(trait) mcmc_new[ , grepl( paste("pred_", trait, "_R_cor.", sep="") , colnames( mcmc_new ) ) ] )

    A<-lapply(1:length(current_vals[[1]][[1]]), function(trait) mcmc_new[ , grepl( paste("pred_", trait, "_A.", sep="") , colnames( mcmc_new ) ) ] )




    ## Define the columns correspondent to the matrix:

    return( list=list(traits=traits, R_sd=R_sd, R_cor=R_cor, A=A) )

  }



  #####tests########

  BM_test<-function(TruePars_scale, phylo, dist="norm", n, save.to.file=T, filename=paste("BM_test.pdf"), quantile=.75){
    #added prior to likelihoodcalc

    if (save.to.file==T){
      pdf(file= paste(filename))
      par(mfrow = c(4, 2))

    }



    lnl_sdw1<-list()
    lnl_sdw2<-list()
    lnl_sdc1<-list()
    lnl_sdc2<-list()
    lnl_Ac2<-list()
    lnl_Ac1<-list()
    lnl_Aw2<-list()
    lnl_Aw1<-list()


    for (w in 1:length(TruePars_scale)){
      print(w)

      tree<-phylo[[w]]

      #widths
      {


        #full_true_res<-lapply(1:2, function(pred) TruePars_scale[[w]]$sim_dat$sim_td_bt[[pred]] )
        #tree<-phylo
        ######R SD w 1#####


        { #R SD w 1

          pred_1_res <- TruePars_scale[[w]]$sim_dat$sim_td[[1]]
          pred_1_A   <- TruePars_scale[[w]]$A$A_ft[[1]]
          pred_1_Rsd <- TruePars_scale[[w]]$R$R_sd[[1]]
          pred_1_Rcor <-TruePars_scale[[w]]$R$R_cor[[1]]

          if(dist=="unif"){

            width_range<-seq(Prior_scale[[1]]$pars$par.sd[2,1],Prior_scale[[1]]$pars$par.sd[2,2], by=.001)
          }else{
            width_range<-sort(rlnorm(n, Prior_scale[[1]]$pars$par.sd[2,1],Prior_scale[[1]]$pars$par.sd[2,2]))
          }
          pred_1_R <- rebuild.cov( r=pred_1_Rcor, v=pred_1_Rsd^2)


          lnl_sdw1[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_1_Rsd[[2]]<-width_range[[x]];
            pred_1_R <- rebuild.cov( r=pred_1_Rcor, v=pred_1_Rsd^2);
            sum( c(
              unlist(lnL_ratematrix(pred_1_res, tree, pred_1_A, pred_1_R, H_fixed = H_fixed)),
              Prior_scale[[1]]$mean.prior(backTransform1(pred_1_A)[1:2]),
              Prior_scale[[1]]$sd.prior(pred_1_Rsd),
              Prior_scale[[1]]$corr.prior(pred_1_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_sdw1[[w]]>quantile(lnl_sdw1[[w]], quantile)], lnl_sdw1[[w]][lnl_sdw1[[w]]>quantile(lnl_sdw1[[w]], quantile)])


          plot(top_sample, main= paste(w, " BM_r_sd_w_",1))
          abline(v=TruePars_scale[[w]]$R$R_sd[[1]][[2]], col=3, lwd=2)
          abline(v=width_range[lnl_sdw1[[w]]==max(lnl_sdw1[[w]])], col=2)


        }



        ######R SD w 2#####

        {

          pred_2_res<-  TruePars_scale[[w]]$sim_dat$sim_td[[2]]
          pred_2_A <-   TruePars_scale[[w]]$A$A_ft[[2]]
          pred_2_Rsd <- TruePars_scale[[w]]$R$R_sd[[2]]
          pred_2_Rcor <-TruePars_scale[[w]]$R$R_cor[[2]]


          if(dist=="unif"){

            width_range<-seq(Prior_scale[[2]]$pars$par.sd[2,1],Prior_scale[[1]]$pars$par.sd[2,2], by=.001)
          }else{
            width_range<-sort(rlnorm(n, Prior_scale[[2]]$pars$par.sd[2,1],Prior_scale[[1]]$pars$par.sd[2,2]))
          }

          pred_2_R <- rebuild.cov( r=pred_2_Rcor, v=pred_2_Rsd^2)


          lnl_sdw2[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_2_Rsd[[2]]<-width_range[[x]];
            pred_2_R <- rebuild.cov( r=pred_2_Rcor, v=pred_2_Rsd^2);
            sum( c(
              unlist(lnL_ratematrix(pred_2_res, tree, pred_2_A, pred_2_R, H_fixed = H_fixed)),
              Prior_scale[[2]]$mean.prior(backTransform1(pred_2_A)[1:2]),
              Prior_scale[[2]]$sd.prior(pred_2_Rsd),
              Prior_scale[[2]]$corr.prior(pred_2_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_sdw2[[w]]>quantile(lnl_sdw2[[w]], quantile)], lnl_sdw2[[w]][lnl_sdw2[[w]]>quantile(lnl_sdw2[[w]], quantile)])


          plot(top_sample, main= paste("BM_r_sd_w_",2))
          abline(v=TruePars_scale[[w]]$R$R_sd[[2]][[2]], col=3, lwd=2)
          abline(v=width_range[lnl_sdw2[[w]]==max(lnl_sdw2[[w]])], col=2)

        }






        ######A W 1#####


        { # A 1

          pred_1_res<-  TruePars_scale[[w]]$sim_dat$sim_td[[1]]
          pred_1_A <-   TruePars_scale[[w]]$A$A_ft[[1]]
          pred_1_Rsd <- TruePars_scale[[w]]$R$R_sd[[1]]
          pred_1_Rcor <-TruePars_scale[[w]]$R$R_cor[[1]]


          if(dist=="unif"){

            width_range<-seq( log(Prior_scale[[1]]$pars$par.mu[2,1]), log(Prior_scale[[1]]$pars$par.mu[2,2]  ), by=.001)
          }else{
            width_range<-sort(log(rlnorm(n, Prior_scale[[1]]$pars$par.mu[2,1],Prior_scale[[1]]$pars$par.mu[2,2])))
          }

          pred_1_R <- rebuild.cov( r=pred_1_Rcor, v=pred_1_Rsd^2)


          lnl_Aw1[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_1_A[[2]]<-width_range[[x]];
            sum( c(
              unlist(lnL_ratematrix(pred_1_res, tree, pred_1_A, pred_1_R, H_fixed = H_fixed)),
              Prior_scale[[1]]$mean.prior(backTransform1(pred_1_A)[1:2]),
              Prior_scale[[1]]$sd.prior(pred_1_Rsd),
              Prior_scale[[1]]$corr.prior(pred_1_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_Aw1[[w]]>quantile(lnl_Aw1[[w]], quantile)], lnl_Aw1[[w]][lnl_Aw1[[w]]>quantile(lnl_Aw1[[w]], quantile)])


          plot(top_sample, main= paste("BM_A_w_",1) )
          abline(v=TruePars_scale[[w]]$A$A_ft[[1]][[2]], col=3, lwd=2)
          abline(v=width_range[lnl_Aw1[[w]]==max(lnl_Aw1[[w]])], col=2)


        }


        ######A W 2#####

        {

          pred_2_res<-  TruePars_scale[[w]]$sim_dat$sim_td[[2]]
          pred_2_A <-   TruePars_scale[[w]]$A$A_ft[[2]]
          pred_2_Rsd <- TruePars_scale[[w]]$R$R_sd[[2]]
          pred_2_Rcor <-TruePars_scale[[w]]$R$R_cor[[2]]

          if(dist=="unif"){

            width_range<-seq( log(Prior_scale[[2]]$pars$par.mu[2,1]), log(Prior_scale[[2]]$pars$par.mu[2,2]  ), by=.001)
          }else{
            width_range<-sort(log(rlnorm(n, Prior_scale[[2]]$pars$par.mu[2,1],Prior_scale[[2]]$pars$par.mu[2,2])))
          }



          pred_2_R <- rebuild.cov( r=pred_2_Rcor, v=pred_2_Rsd^2)

          lnl_Aw2[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_2_A[[2]]<-width_range[[x]];
            sum( c(
              unlist(lnL_ratematrix(pred_2_res, tree, pred_2_A, pred_2_R, H_fixed = H_fixed)),
              Prior_scale[[2]]$mean.prior(backTransform1(pred_2_A)[1:2]),
              Prior_scale[[2]]$sd.prior(pred_2_Rsd),
              Prior_scale[[2]]$corr.prior(pred_2_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_Aw2[[w]]>quantile(lnl_Aw2[[w]], quantile)], lnl_Aw2[[w]][lnl_Aw2[[w]]>quantile(lnl_Aw2[[w]], quantile)])


          plot(top_sample, main= paste("BM_A_w_",2) )
          abline(v=TruePars_scale[[w]]$A$A_ft[[2]][[2]], col=3, lwd=2)
          abline(v=width_range[lnl_Aw2[[w]]==max(lnl_Aw2[[w]])], col=2)

        }

      }


      #centers
      {

        ######R SD c 1#####


        {

          pred_1_res <- TruePars_scale[[w]]$sim_dat$sim_td[[1]]
          pred_1_A   <- TruePars_scale[[w]]$A$A_ft[[1]]
          pred_1_Rsd <- TruePars_scale[[w]]$R$R_sd[[1]]
          pred_1_Rcor <-TruePars_scale[[w]]$R$R_cor[[1]]


          if(dist=="unif"){

            width_range<-seq(Prior_scale[[1]]$pars$par.sd[1,1],Prior_scale[[1]]$pars$par.sd[1,2], by=.001)


          }else{
            width_range<-sort(rlnorm(n, Prior_scale[[1]]$pars$par.sd[1,1],Prior_scale[[1]]$pars$par.sd[1,2]))
          }



          pred_1_R <- rebuild.cov( r=pred_1_Rcor, v=pred_1_Rsd^2)




          lnl_sdc1[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_1_Rsd[[1]]<-width_range[[x]];
            pred_1_R <- rebuild.cov( r=pred_1_Rcor, v=pred_1_Rsd^2);
            sum( c(
              unlist(lnL_ratematrix(pred_1_res, tree, pred_1_A, pred_1_R, H_fixed = H_fixed)),
              Prior_scale[[1]]$mean.prior(backTransform1(pred_1_A)[1:2]),
              Prior_scale[[1]]$sd.prior(pred_1_Rsd),
              Prior_scale[[1]]$corr.prior(pred_1_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_sdc1[[w]]>quantile(lnl_sdc1[[w]], quantile)], lnl_sdc1[[w]][lnl_sdc1[[w]]>quantile(lnl_sdc1[[w]], quantile)])


          plot(top_sample, main= paste(w, " BM_r_sd_c_",1))
          abline(v=TruePars_scale[[w]]$R$R_sd[[1]][[1]], col=3, lwd=2)
          abline(v=width_range[lnl_sdc1[[w]]==max(lnl_sdc1[[w]])], col=2)


        }



        ######R SD c 2#####

        {


          pred_2_res<-  TruePars_scale[[w]]$sim_dat$sim_td[[2]]
          pred_2_A <-   TruePars_scale[[w]]$A$A_ft[[2]]
          pred_2_Rsd <- TruePars_scale[[w]]$R$R_sd[[2]]
          pred_2_Rcor <-TruePars_scale[[w]]$R$R_cor[[2]]

          if(dist=="unif"){


            width_range<-seq(Prior_scale[[2]]$pars$par.sd[1,1],Prior_scale[[2]]$pars$par.sd[1,2], by=.001)

          }else{
            width_range<-sort(rlnorm(n, Prior_scale[[2]]$pars$par.sd[1,1],Prior_scale[[2]]$pars$par.sd[1,2]))
          }



          pred_2_R <- rebuild.cov( r=pred_2_Rcor, v=pred_2_Rsd^2)


          lnl_sdc2[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_2_Rsd[[1]]<-width_range[[x]];
            pred_2_R <- rebuild.cov( r=pred_2_Rcor, v=pred_2_Rsd^2);
            sum( c(
              unlist(lnL_ratematrix(pred_2_res, tree, pred_2_A, pred_2_R, H_fixed = H_fixed)),
              Prior_scale[[2]]$mean.prior(backTransform1(pred_2_A)[1:2]),
              Prior_scale[[2]]$sd.prior(pred_2_Rsd),
              Prior_scale[[2]]$corr.prior(pred_2_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_sdc2[[w]]>quantile(lnl_sdc2[[w]], quantile)], lnl_sdc2[[w]][lnl_sdc2[[w]]>quantile(lnl_sdc2[[w]], quantile)])


          plot(top_sample, main= paste("BM_r_sd_c_",2))
          abline(v=TruePars_scale[[w]]$R$R_sd[[2]][[1]], col=3, lwd=2)
          abline(v=width_range[lnl_sdc2[[w]]==max(lnl_sdc2[[w]])], col=2)

        }




        ######A C 1#####


        { # A 1

          pred_1_res<-  TruePars_scale[[w]]$sim_dat$sim_td[[1]]
          pred_1_A <-   TruePars_scale[[w]]$A$A_ft[[1]]
          pred_1_Rsd <- TruePars_scale[[w]]$R$R_sd[[1]]
          pred_1_Rcor <-TruePars_scale[[w]]$R$R_cor[[1]]



          if(dist=="unif"){


            width_range<-seq( Prior_scale[[1]]$pars$par.mu[1,1] ,Prior_scale[[1]]$pars$par.mu[1,2], by=.001)

          }else{
            width_range<-sort(rnorm(n, Prior_scale[[1]]$pars$par.mu[1,1],Prior_scale[[1]]$pars$par.mu[1,2]))
          }







          pred_1_R <- rebuild.cov( r=pred_1_Rcor, v=pred_1_Rsd^2)


          lnl_Ac1[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_1_A[[1]]<-width_range[[x]];
            sum( c(
              unlist(lnL_ratematrix(pred_1_res, tree, pred_1_A, pred_1_R, H_fixed = H_fixed)),
              Prior_scale[[1]]$mean.prior(backTransform1(pred_1_A)[1:2]),
              Prior_scale[[1]]$sd.prior(pred_1_Rsd),
              Prior_scale[[1]]$corr.prior(pred_1_Rcor)))
          }
          ))

          top_sample<-cbind(width_range[lnl_Ac1[[w]]>quantile(lnl_Ac1[[w]], quantile)], lnl_Ac1[[w]][lnl_Ac1[[w]]>quantile(lnl_Ac1[[w]], quantile)])


          plot(top_sample, main= paste("BM_A_c_",1))
          abline(v=TruePars_scale[[w]]$A$A_ft[[1]][[1]], col=3, lwd=2)
          abline(v=width_range[lnl_Ac1[[w]]==max(lnl_Ac1[[w]])], col=2)

          #phenogram(phylo[[w]], setNames(unlist(pred_1_res[,1]), phylo[[w]]$tip.label), ftype = "off"  )
          #abline(h=TruePars_scale[[w]]$A$A_ft[[1]][1], col=3)
          #abline(h=width_range[lnl_Ac1[[w]]==max(lnl_Ac1[[w]])], col=2)


        }


        ######A C 2#####

        {

          pred_2_res<-  TruePars_scale[[w]]$sim_dat$sim_td[[2]]
          pred_2_A <-   TruePars_scale[[w]]$A$A_ft[[2]]
          pred_2_Rsd <- TruePars_scale[[w]]$R$R_sd[[2]]
          pred_2_Rcor <-TruePars_scale[[w]]$R$R_cor[[2]]




          if(dist=="unif"){


            width_range<-seq( Prior_scale[[2]]$pars$par.mu[1,1] ,Prior_scale[[2]]$pars$par.mu[1,2], by=.001)

          }else{
            width_range<-sort(rnorm(n, Prior_scale[[2]]$pars$par.mu[1,1],Prior_scale[[2]]$pars$par.mu[1,2]))
          }


          pred_2_R <- rebuild.cov( r=pred_2_Rcor, v=pred_2_Rsd^2)
          #TruePars_scale[[1]]$R$R

          lnl_Ac2[[w]]<-unlist(lapply(1:length(width_range), function(x) {
            pred_2_A[[1]]<-width_range[[x]];
            sum( c(
              unlist(lnL_ratematrix(pred_2_res, tree, pred_2_A, pred_2_R, H_fixed = H_fixed)),
              Prior_scale[[2]]$mean.prior(backTransform1(pred_2_A)[1:2]),
              Prior_scale[[2]]$sd.prior(pred_2_Rsd),
              Prior_scale[[2]]$corr.prior(pred_2_Rcor)))
          }
          ))


          top_sample<-cbind(width_range[lnl_Ac2[[w]]>quantile(lnl_Ac2[[w]], quantile)], lnl_Ac2[[w]][lnl_Ac2[[w]]>quantile(lnl_Ac2[[w]], quantile)])



          plot(top_sample, main= paste("BM_A_c_",2))
          abline(v=TruePars_scale[[w]]$A$A_ft[[2]][[1]], col=3, lwd=2)
          abline(v=width_range[lnl_Ac2[[w]]==max(lnl_Ac2[[w]])], col=2)

        }

      }

    }
    if (save.to.file==T){
      dev.off()
    }

  }



}
