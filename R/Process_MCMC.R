

read_BePhyNE_log=function(file_name){
  read.table(file_name,sep = "\t", header = T)
}

colMedians_df <- function(df) {
  is_num <- sapply(df, is.numeric)
  medians <- sapply(df[is_num], median, na.rm = TRUE)
  as.data.frame(t(medians), stringsAsFactors = FALSE)
}


logdf2traitsdf_list = function(logdf, transform2nichespace=T, logrows= 1:nrow(logdf)){

  predictors <- unique(gsub("pred_(\\d+)_.*", "\\1", grep("^pred_\\d+_", names(logdf), value = TRUE)))
  pred_prefix <- paste0("pred_", predictors, "_")

  traits_df_list= list() #lapply(predictors, function(pred) list())


  count=0
  for(x in logrows){
  count=count+1
    for(pred in predictors){

      pred_prefix <- paste0("pred_", pred, "_")

      print(x)
      # --- Trait matrix ---
      trait_cols <- grep(paste0("^", pred_prefix, "dat\\.trait"), names(logdf), value = TRUE)
      trait_mat <- logdf[x, trait_cols, drop = FALSE]
      #trait_medians <- apply(trait_mat, 2, median, na.rm = TRUE)

      # Extract row and column from names like "pred_1_dat.traitXYZ"
      trait_df <- do.call(rbind, lapply(names(trait_mat), function(name) {
        num <- as.integer(gsub(".*trait", "", name))
        if ( num<100){
          col <- num %/% 10
          row <- num %% 10
        }else{
          col <- num %/% 100
          row <- num %% 100
        }
        if (row == 0) {
          col <- col - 1
          row <- 100
        }
        data.frame(row = row, col = col, value = trait_mat[[name]])
      }))

      # Create matrix: rows = 1:82, cols = 1:ncol
      nrow_out <- max(trait_df$row)
      ncol_out <- max(trait_df$col)
      trait_matrix <- matrix(NA_real_, nrow = nrow_out, ncol = ncol_out)
      for (i in seq_len(nrow(trait_df))) {
        trait_matrix[trait_df$row[i], trait_df$col[i]] <- trait_df$value[i]
      }
      rownames(trait_matrix) <- 1:nrow_out
      colnames(trait_matrix) <- c(paste0("opt"), paste0("brdth"), paste0("tol"))

      if(transform2nichespace==T){
        BT_trait_matrix = do.call(rbind, lapply(1:nrow(trait_matrix), function(i) backTransform1(trait_matrix[i,])))
        traits_df_list[[pred]][[count]] = BT_trait_matrix
      }else{
        traits_df_list[[pred]][[count]] = trait_matrix

      }

    }
  }
  return(traits_df_list)
}



logdf2R_list = function(logdf, logrows= 1:nrow(logdf)){


  predictors <- unique(gsub("pred_(\\d+)_.*", "\\1", grep("^pred_\\d+_", names(logdf), value = TRUE)))


  Rsd_list  = list() #lapply(predictors, function(pred) list())
  Rcor_list = list() #lapply(predictors, function(pred) list())
  R_list    = list() #lapply(predictors, function(pred) list())

  count=0
  
  for(i in logrows){
    count=count+1
    
    for(pred in predictors){

      pred_prefix <- paste0("pred_", pred, "_")

      # --- R_sd vector ---
      Rsd_cols <- grep(paste0("^", pred_prefix, "R_sd\\d+"), names(logdf), value = TRUE)
      Rsd <- unlist(logdf[i,Rsd_cols])
      #Rsd_median <- unlist(apply(Rsd_mat, 2, median, na.rm = TRUE))

      # --- R_cor vector ---
      Rcor_cols <- grep(paste0("^", pred_prefix, "R_cor\\d+"), names(logdf), value = TRUE)
      Rcor      <-  matrix(as.numeric(logdf[i,Rcor_cols]), nrow = length(Rsd), byrow = TRUE)
      #Rcor_median <- matrix(apply(Rcor_mat, 2, median, na.rm = TRUE), nrow = length(Rsd_median), byrow = TRUE)

      # --- Rate matrix ---
      #Rval_cols <- grep(paste0("^", pred_prefix, "R\\d+"), names(logdf), value = TRUE)
      #Rval_mat <- logdf[ Rval_cols]
      #R_median <- matrix(apply(Rval_mat, 2, median, na.rm = TRUE), nrow = length(Rsd_median), byrow = TRUE)

      R = corpcor::rebuild.cov( r = Rcor, v = as.numeric(as.vector(Rsd)))

       Rsd_list[[pred]][[count]]  = Rsd
      Rcor_list[[pred]][[count]] = Rcor
         R_list[[pred]][[count]]    = R

    }
  }

  return(list( Rsd = Rsd_list
              ,Rcor= Rcor_list
              ,R   = R_list
              )
         )

}


logdf2A_list = function(logdf, transform2nichespace=T, logrows= 1:nrow(logdf)){

  predictors <- unique(gsub("pred_(\\d+)_.*", "\\1", grep("^pred_\\d+_", names(logdf), value = TRUE)))

  A_list  =  list() #lapply(predictors, function(pred) list())

  count=0
  for(i in logrows){
    count=count+1
    for(pred in predictors){

      pred_prefix <- paste0("pred_", pred, "_")

      A_cols <- grep(paste0("^", pred_prefix, "A\\d+"), names(logdf), value = TRUE)
      A_mat <- logdf[i,A_cols]

      if(transform2nichespace==T){
        A_list[[pred]][[count]]  = backTransform1(unlist(c(A_mat,1) ))[1:2]
      }else{

        A_list[[pred]][[count]]  = A_mat
      }
    }
  }

  return( A_list
  )

}


logdf2parlist = function(logdf, transform2nichespace=T, logrows= 1:nrow(logdf)){

  traits_list = logdf2traitsdf_list(logdf, transform2nichespace, logrows)

  A_list  = logdf2A_list(logdf, transform2nichespace, logrows)

  R_lists = logdf2R_list(logdf, logrows)
  R_lists$Rsd
  R_lists$Rcor
  R_lists$R

  return(list(traits = traits_list
              ,A     = A_list
              ,Rsd   = R_lists$Rsd
              ,Rcor  = R_lists$Rcor
              ,R     = R_lists$R

  ))
}


logdf2medians = function(logdf){

  medians = colMedians_df(logdf)

  return(logdf2parlist(medians))


}



summarize_logdf=function(logdf, HPD_prob=0.95){
  
  mcmc_obj = coda::as.mcmc(logdf[,-1])
  mcmc_ESS = effectiveSize(mcmc_obj)
  mcmc_HPD = HPDinterval(mcmc_obj, prob = HPD_prob)

  mcmc_median_df      = colMedians_df(logdf)
  mcmc_median_parlist = logdf2medians(logdf)

  summary = list(ESS = mcmc_ESS
                 ,HPD = mcmc_HPD
                 ,median_df = mcmc_median_df
                 ,median_parlist = mcmc_median_parlist
                 )
}
