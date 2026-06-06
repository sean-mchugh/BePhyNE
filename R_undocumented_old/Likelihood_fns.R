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
    #cur.liks$glm.lik[missSP]
    prop.glm.lik[missSP] <- missing.glm.lik
    prop.glm.sum.lik <- sum(unlist(prop.glm.lik))
    GLM <- prop.glm.sum.lik-cur.liks$glm.sum.lik


    if(prop==0){

      #if height was chosen as the prior then you must calculate how the prior changes



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



