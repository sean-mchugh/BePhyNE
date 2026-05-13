



make_simmaps_BePhyNE=function( tree, 
                               logdf,
                               char_names =NA, 
                               nsims = 10){
  
  parlist = logdf2parlist(logdf, transform2nichespace = F, logrows = sample(1:nrow(logdf),nsims))
  
  if(sum(is.na(char_names))==1){
    
     char_names = unlist(lapply(1:length(parlist$A), function(i) c(paste0("X",i ,"_optimum"), paste0("X",i ,"_breadth"))))
    
  }

  traits_list = lapply(1:length(parlist$traits[[1]]), function(i) do.call(cbind, lapply(parlist$traits, function(pred) pred[[i]][,1:2])))
  
  R_list = lapply(1:length(parlist$R[[1]]), function(i) as.matrix(do.call(bdiag, lapply(parlist$R, function(pred) pred[[i]]))))
  

  for (i in 1:length(traits_list)){
    
    colnames(traits_list[[i]])  = char_names
    rownames(traits_list[[i]])  = tree$tip.label
    colnames(R_list[[i]])       = char_names
  }
  
  sims = contsimmap::make.contsimmap(tree,  trait.data = traits_list, Xsig2 = R_list,  nsims = nsims)
  
  return(sims)
}
