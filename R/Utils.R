
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

make_tuning <- function(tree, pred,
                        center_slide = 0.18,
                        center_mult  = 0.12,
                        width_slide  = 0.15,
                        width_mult   = 0.2,
                        height_slide = 0.5,
                        height_mult  = 0.3,
                        w_mu  = c(0.6, 0.6),
                        w_sd  = c(0.15, 0.12),
                        v_cor = 100) {
  
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
  
  return(tuning)
}