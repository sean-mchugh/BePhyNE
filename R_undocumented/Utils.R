
###Utils##################################################################################################################################################################

'%!in%' <- function(x,y)!('%in%'(x,y))

format_BePhyNE_data = function(pa_data, tree, sp_col, occ_col, env_preds, scale_atr=NA){

  pres_data_scaled= pa_data

  if( sum(is.na(scale_atr))>0){
    scaled_climate = scale(pa_data[,env_preds])

    pres_data_scaled[,env_preds]= scaled_climate[,env_preds]

    scale_atr <- list(scale=lapply(env_preds, function(pred)(attr(scaled_climate , "scaled:scale")[colnames(scaled_climate)==pred])),
                      center= lapply(env_preds, function(pred) (attr(scaled_climate,  "scaled:center")[colnames(scaled_climate)==pred]))
    )

    names(scale_atr$scale) = env_preds
    names(scale_atr$center) = env_preds


  } else{


    scaled_climate = do.call(cbind ,lapply(env_preds, function(pred) (pa_data[,pred]-scale_atr$center[[pred]])/scale_atr$scale[[pred]]) )

    colnames(scaled_climate) = env_preds

    pres_data_scaled[,env_preds]= scaled_climate[,env_preds]

  }

  data_final <- lapply(split(pres_data_scaled[c(occ_col, env_preds)],pres_data_scaled[sp_col]), as.list)
  #replace vector of species names with a single name
  data_final_tree = list()
  for (sp in 1:length(tree$tip.label)) {
    sp_name= tree$tip.label[[sp]]

    data_final[[sp_name]]$species= names(data_final[sp_name])[[1]]
    data_final[[sp_name]] = data_final[[sp_name]][c(sp_col, occ_col, env_preds)]
    names( data_final[[sp_name]]) =c("species", "y", paste0("X", 1:(length(names( data_final[[sp_name]]))-2)) )

    data_final_tree[[sp_name]] = data_final[[sp_name]]
  }

  # match data_final species order to tree tip order
  #data_final=data_final[order(match(unlist(lapply(data_final, function(sp) sp[[sp_col]])), tree$tip.label))]

  return(list(data=data_final_tree, scale = scale_atr))
}


###Utils##################################################################################################################################################################

data_df2list=function(df, tree)
{

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

make_tuning_par_obj=function( tree, npred, center_slide=.18,
                              center_mult=.12,
                              width_slide=.15,
                              width_mult=.2,
                              height_slide=.5,
                              height_mult=.3,
                              w_mu_slide=c(0.6,0.6),
                              w_mu_multi=c(0.6,0.6),
                              w_sd_slide=c(.15,.12),
                              w_sd_multi=c(.15,.12),
                              v_cor_iwish=100)
{
  pred=npred
  tuning<-  list(
    niche_prop= lapply(1:pred, function(pred) list(slide=tibble(center=sample(center_slide, length(tree$tip.label), replace=T),
                                                                width= sample(width_slide, length(tree$tip.label), replace=T),
                                                                height= sample(height_slide, length(tree$tip.label), replace=T) ),
                                                   multi=tibble(center=sample(center_mult , length(tree$tip.label), replace=T),
                                                                width= sample(width_mult , length(tree$tip.label), replace=T),
                                                                height= sample(height_mult, length(tree$tip.label), replace=T)  )
    )),
    w_mu =  lapply(1:pred, function(pred) list(slide=w_mu_slide,
                                               multi=w_mu_multi))
    ,
    w_sd =  lapply(1:pred, function(pred)  list(slide=w_sd_slide,
                                                multi=w_sd_multi)
    ),
    v_cor       = lapply(1:pred, function(pred) v_cor_iwish)
  )

  return(tuning)

}
