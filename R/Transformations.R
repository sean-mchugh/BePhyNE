##transformations##################################################################################################################################################################

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




backtransform_denormalize_parlist = function(parlist,
                                             scale_atr,
                                             only_denormalize = TRUE) {
  
  out = parlist
  
  for (pred in names(out$traits)) {
    
    pred_int = as.integer(pred)
    pred_scale  = unname(scale_atr$scale[pred_int])
    pred_center = unname(scale_atr$center[pred_int])
    
    out$traits[[pred]] = lapply(out$traits[[pred]], function(mat) {
      
      if (!only_denormalize) {
        mat[, "brdth"] = exp(mat[, "brdth"])
        mat[, "tol"] = 0.05 + 0.95 * (exp(-mat[, "tol"]) / (1 + exp(-mat[, "tol"])))
      }
      
      mat[, "opt"]   = mat[, "opt"] * pred_scale + pred_center
      mat[, "brdth"] = mat[, "brdth"] * pred_scale 
      
      mat
    })
    
    out$A[[pred]] = lapply(out$A[[pred]], function(a) {
      
      A1 = paste0("pred_", pred, "_A1")
      A2 = paste0("pred_", pred, "_A2")
      
      if (!only_denormalize) {
        a[A2] = exp(a[A2])
      }
      
      a[A1] = a[A1] * pred_scale + pred_center
      a[A2] = a[A2] * pred_scale 
      
      a
    })
    
    out$Rsd[[pred]] = lapply(out$Rsd[[pred]], function(R_sd) {
      R_sd * pred_scale
    })
    
    out$R[[pred]] = Map(
      function(R_sd, prop_cor) {
        corpcor::rebuild.cov(r = prop_cor, v = R_sd^2)
      },
      out$Rsd[[pred]],
      out$Rcor[[pred]]
    )
  }
  
  out
}

backtransform_denormalize_parlist = function(parlist,
                                             scale_atr = NA,
                                             only_denormalize = TRUE) {
  
  out = parlist
  do_denorm = !all(is.na(scale_atr))
  
  for (pred in names(out$traits)) {
    
    if (do_denorm) {
      pred_int = as.integer(pred)
      pred_scale  = unname(scale_atr$scale[pred_int])
      pred_center = unname(scale_atr$center[pred_int])
    } else {
      pred_scale  = 1
      pred_center = 0
    }
    
    out$traits[[pred]] = lapply(out$traits[[pred]], function(mat) {
      
      if (!only_denormalize) {
        mat[, "brdth"] = exp(mat[, "brdth"])
        mat[, "tol"] = 0.05 + 0.95 * (
          exp(-mat[, "tol"]) / (1 + exp(-mat[, "tol"]))
        )
      }
      
      if (do_denorm) {
        mat[, "opt"]   = mat[, "opt"] * pred_scale + pred_center
        mat[, "brdth"] = mat[, "brdth"] * pred_scale 
      }
      
      ## tolerance never denormalized
      
      mat
    })
    
    out$A[[pred]] = lapply(out$A[[pred]], function(a) {
      
      A1 = paste0("pred_", pred, "_A1")
      A2 = paste0("pred_", pred, "_A2")
      
      if (!only_denormalize) {
        a[A2] = exp(a[A2])
      }
      
      if (do_denorm) {
        a[A1] = a[A1] * pred_scale + pred_center
        a[A2] = a[A2] * pred_scale 
      }
      
      a
    })
    
    if (do_denorm) {
      
      out$Rsd[[pred]] = lapply(out$Rsd[[pred]], function(R_sd) {
        R_sd * pred_scale
      })
      
      out$R[[pred]] = Map(
        function(R_sd, prop_cor) {
          corpcor::rebuild.cov(
            r = prop_cor,
            v = R_sd^2
          )
        },
        out$Rsd[[pred]],
        out$Rcor[[pred]]
      )
    }
  }
  
  out
}

backtransform_denormalize_logdf = function(logdf, scale_atr = NA) {
  
  out = logdf
  do_denorm = !all(is.na(scale_atr))
  
  pred_ids = as.integer(sub(
    "^pred_(\\d+).*$", "\\1",
    grep("^pred_\\d+_", names(out), value = TRUE)
  ))
  pred_ids = sort(unique(pred_ids))
  
  for (pred in pred_ids) {
    
    if (do_denorm) {
      pred_scale  = unname(scale_atr$scale[pred])
      pred_center = unname(scale_atr$center[pred])
    } else {
      pred_scale  = 1
      pred_center = 0
    }
    
    opt_cols   = grep(paste0("^pred_", pred, "_dat\\.opt_"), names(out), value = TRUE)
    brdth_cols = grep(paste0("^pred_", pred, "_dat\\.brdth_"), names(out), value = TRUE)
    tol_cols   = grep(paste0("^pred_", pred, "_dat\\.tol_"), names(out), value = TRUE)
    
    if (length(opt_cols) > 0) {
      out[opt_cols] = out[opt_cols] * pred_scale + pred_center
    }
    
    if (length(brdth_cols) > 0) {
      out[brdth_cols] = exp(as.matrix(out[brdth_cols])) * pred_scale 
    }
    
    if (length(tol_cols) > 0) {
      x = as.matrix(out[tol_cols])
      out[tol_cols] = 0.05 + 0.95 * (exp(-x) / (1 + exp(-x)))
    }
    
    A1 = paste0("pred_", pred, "_A1")
    A2 = paste0("pred_", pred, "_A2")
    
    if (A1 %in% names(out)) {
      out[[A1]] = out[[A1]] * pred_scale + pred_center
    }
    
    if (A2 %in% names(out)) {
      out[[A2]] = exp(out[[A2]]) * pred_scale
    }
    
    if (do_denorm) {
      
      Rsd1 = paste0("pred_", pred, "_R_sd1")
      Rsd2 = paste0("pred_", pred, "_R_sd2")
      
      out[[Rsd1]] = out[[Rsd1]] * pred_scale
      out[[Rsd2]] = out[[Rsd2]] * pred_scale
      
      R_cols = paste0("pred_", pred, "_R", 1:4)
      Rcor_cols = paste0("pred_", pred, "_R_cor", 1:4)
      
      for (i in seq_len(nrow(out))) {
        
        R_sd = c(out[[Rsd1]][i], out[[Rsd2]][i])
        
        prop_cor = matrix(
          as.numeric(out[i, Rcor_cols]),
          nrow = 2,
          ncol = 2,
          byrow = TRUE
        )
        
        out[i, R_cols] = as.numeric(corpcor::rebuild.cov(
          r = prop_cor,
          v = R_sd^2
        ))
      }
    }
  }
  
  out
}


backtransform_denormalize_logsummary = function(log_summary, scale_atr=NA) {
  
  out = log_summary
  
  out$median_df = backtransform_denormalize_logdf(
    logdf = out$median_df,
    scale_atr = scale_atr
  )
  
  out$median_parlist = backtransform_denormalize_parlist(
    parlist = out$median_parlist,
    scale_atr = scale_atr,
    only_denormalize = F
  )
  
  out$HPDlower_parlist = backtransform_denormalize_parlist(
    parlist = out$HPDlower_parlist,
    scale_atr = scale_atr,
    only_denormalize = F
  )
  
  out$HPDupper_parlist = backtransform_denormalize_parlist(
    parlist = out$HPDupper_parlist,
    scale_atr = scale_atr,
    only_denormalize = F
  )
  
  out
}


