
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


plot_response_curves <- function(tree,
                                 models_list,       # list of models -> each is a list of matrices: one per predictor
                                 model_names    = paste0("model", 1:length(models_list)),
                                 predictor_names= paste0("pred", 1:length(models_list[[1]])),   # character vector of predictor names
                                 scale_atr=NA,         # list with $center and $scale (same length as predictors)
                                 curve_colors      = rep(make.transparent("black", 255/255), length(models_list)),
                                 curve_fill_colors = rep(make.transparent("white", 255/255), length(models_list)),
                                 line_types    = rep(1, length(models_list)),
                                 xlims = NULL
) {
  
  
  #if(class(models_list[[1]][[1]])!="list"){
  #  models_list = list(models_list)
  #}
  
  
  library(ape)
  response_function = function(x) 1 / (1 + exp(-x))
  
  npreds <- length(predictor_names)
  nmodels <- length(models_list)
  ntips <- length(tree$tip.label)
  
  #if (is.null(curve_colors)) curve_colors <- rep("black", nmodels)
  #if (is.null(line_types)) line_types <- rep(1, nmodels)
  
  tree$edge.length <- tree$edge.length / max(branching.times(tree))
  TL <- max(branching.times(tree))
  xplots <- lapply(1:npreds, function(i) c((2*i)*TL, (2*i + 1.5)*TL))
  maxX <- max(sapply(xplots, max)) + 0.3*TL
  
  # Compute predictor sequences from θ ± 3ω for all models
  predictor_ranges <- vector("list", npreds)
  predictor_sequences <- vector("list", npreds)
  
  #for (p in seq_len(npreds)) {
  #  all_opt <- c()
  #  all_breadth <- c()
  #  
  #  for (model in models_list) {
  #    mat <- model[[p]]  # θ, ω, t
  #    all_opt <- c(all_opt, mat[,1])
  #    all_breadth <- c(all_breadth, mat[,2])
  #  }
  #  lower <- min(all_opt - all_breadth)
  #  upper <- max(all_opt + all_breadth)
  #  predictor_ranges[[p]] <- c(lower, upper)
  #  predictor_sequences[[p]] <- seq(lower, upper, length.out = 200)
  #}
  
  
  predictor_ranges <- vector("list", npreds)
  predictor_sequences <- vector("list", npreds)
  
  for (p in seq_len(npreds)) {
    if (!is.null(xlims) && !is.null(xlims[[p]])) {
      predictor_ranges[[p]] <- xlims[[p]]
    } else {
      all_opt <- c()
      all_breadth <- c()
      for (model in models_list) {
        mat <- model[[p]]  # θ, ω, t
        all_opt <- c(all_opt, mat[,1])
        all_breadth <- c(all_breadth, mat[,2])
      }
      predictor_ranges[[p]] <- range(all_opt - all_breadth, all_opt + all_breadth)
    }
    predictor_sequences[[p]] <- seq(predictor_ranges[[p]][1], predictor_ranges[[p]][2], length.out = 200)
  }
  
  # Convert all model parameters to coefficient matrices
  model_coefs <- lapply(models_list, function(model) {
    lapply(model, function(pred) traits2coefs(pred) )
  })
  
  par(mar = c(0, 1.5, 0, 1.5), mgp = c(0, 0, 0))
  plot.phylo(tree, x.lim = c(0, maxX), y.lim = c(-3, ntips + 3), cex = 0.3, no.margin = TRUE)
  
  for (pred_idx in seq_len(npreds)) {
    xseq <- predictor_sequences[[pred_idx]]
    xplot_range <- xplots[[pred_idx]]
    xxp <- seq(xplot_range[1], xplot_range[2], length.out = length(xseq))
    
    for (tip in 1:ntips) {
      lines(xxp, rep(tip, length(xxp)), col = if (tip %% 2 == 0) "grey" else "grey40", lwd = 1)
      
    }
    
    for (model_idx in seq_along(model_coefs)) {
      coefs <- model_coefs[[model_idx]][[pred_idx]]
      for (tip in 1:ntips) {
        β <- coefs[tip, ]
        yvals <- response_function(β[[1]] + β[[2]]*xseq + β[[3]]*xseq^2)
        
        if(!is.na(curve_fill_colors[[model_idx]])){
          polygon(xxp[yvals >.01], yvals[yvals > 0.01] + tip, col=curve_fill_colors[[model_idx]],density=NA,lwd=.5,border=curve_colors[[model_idx]])
        }
        if(!is.na(line_types[[model_idx]])){
          lines(xxp[yvals > 0.01], yvals[yvals > 0.01] + tip,
                col = curve_colors[model_idx],
                lwd = 0.5,
                lty = line_types[model_idx])
        }
        peak <- which.max(yvals)
        points(xxp[peak], tip, pch = "*", col = curve_colors[model_idx], cex = 0.4)
      }
    }
    
    if(sum(is.na(scale_atr))==1){
      xlabels <- round(seq(predictor_ranges[[pred_idx]][1], predictor_ranges[[pred_idx]][2], length.out = 4), 2)
      
    }else{
      xlabels <- round(seq(predictor_ranges[[pred_idx]][1], predictor_ranges[[pred_idx]][2], length.out = 4) *
                         scale_atr$scale[[pred_idx]] + scale_atr$center[[pred_idx]], 2)
    }
    label_xloc <- seq(xplot_range[1], xplot_range[2], length.out = 4)
    
    text(x = label_xloc, y = rep(-0.2, 4), labels = xlabels, cex = 0.4)
    text(x = mean(xplot_range), y = ntips + 2, labels = predictor_names[pred_idx], cex = 0.4)
  }
  
  rep(NA, length(model_names))
  pch =rep(NA, length(models_list))
  pch[!is.na(curve_fill_colors)] =22
  
  pt.bg =curve_colors
  pt.bg[!is.na(curve_fill_colors)] =curve_fill_colors[!is.na(curve_fill_colors)]
  
  legend(x="bottomright",legend=model_names,pch=pch, lty=line_types, pt.cex=1.5,pt.bg=pt.bg,
         box.col="transparent", cex=.8,ncol=3,)
  
}



plot_summary_ridgeplot = function(tree,
                                  log_summary  ,  
                                  model_names       = paste0("model", 1),
                                  predictor_names   = NA,   # character vector of predictor names
                                  scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                                  curve_colors      = rep(make.transparent("black", 255/255),1),
                                  curve_fill_colors = rep(make.transparent("white", 255/255), 1),
                                  line_types        = rep(1, 1),
                                  xlims = NULL
){
  
  
  models_list = list(lapply(log_summary$median_parlist$traits, function(i) i[[1]]))
  
  if(sum(is.na(predictor_names))==1){
    
    predictor_names = paste0("Predictor ", 1:length(model_list[[1]]))
    
  }
  
  plot_response_curves(tree
                       ,models_list       = models_list       # list of models -> each is a list of matrices: one per predictor
                       ,model_names       = model_names      
                       ,predictor_names   = predictor_names  
                       ,scale_atr         = scale_atr               
                       ,curve_colors      = curve_colors     
                       ,curve_fill_colors = curve_fill_colors
                       ,line_types        = line_types       
                       ,xlims             = xlims 
  )
  
  
  
  
}



plot_summarylist_ridgeplot = function(tree,
                                  log_summarylist  ,  
                                  model_names       = paste0("model", length(log_summarylist)),
                                  predictor_names   = NA,   # character vector of predictor names
                                  scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                                  curve_colors      = rep(make.transparent("black", 255/255), length(log_summarylist)),
                                  curve_fill_colors = rep(make.transparent("white", 255/255), length(log_summarylist)),
                                  line_types        = rep(1, length(log_summarylist)),
                                  xlims = NULL
){
  
  
  model_list = lapply(log_summarylist, function (x) lapply(x$median_parlist$traits, function(i) i[[1]]))
  
  if(sum(is.na(predictor_names))==1){
    
    predictor_names = paste0("Predictor ", 1:length(model_list[[1]]))
    
  }
  
  plot_response_curves(tree
                       ,models_list       = model_list       # list of models -> each is a list of matrices: one per predictor
                       ,model_names       = model_names      
                       ,predictor_names   = predictor_names  
                       ,scale_atr         = scale_atr               
                       ,curve_colors      = curve_colors     
                       ,curve_fill_colors = curve_fill_colors
                       ,line_types        = line_types       
                       ,xlims             = xlims 
  )
  
  
  
  
}


plot_AUC_treebarplot = function(tree, predict_stats_list, cols= c("blue"), xlim =c(0,100), fsize = 0.6, mar = c(5.1, 1, 1.1, 0.5),  label.offset=1){
  
  AUC = unlist(lapply(predict_stats_list, function(i) i$AUC))
  
  names(AUC) = names(predict_stats_list)
  
  
  plotTree.barplot(
    tree,
    AUC,
    args.barplot = list(
      col = cols,
      border = cols,
      xlab = "AUC",
      xlim=xlim,
      mar = c(5.1, 0, 0.1, 4)  # right margin for barplot
    ),
    args.plotTree = list(
      fsize = fsize,                # smaller tip labels
      mar = mar,  # left margin for tree
      label.offset          # more offset between tree and bars
    )
  )
  
}






legend.evorates_mod =function (sim, location = c("bottomleft", "topleft", "bottomright", 
                                                 "topright"), color.element = "R", exp = FALSE, exp.txt = TRUE, 
                               col = def.color.scheme(), val.range = NULL, res = 100, alpha = NA, 
                               breaks = NULL, select.levels = NULL, box.dims = NULL, box.offset = NULL, 
                               box.scale = 1, txt.col = NULL, main = NULL, labels_n=3, ...) 
{
  if (exp) {
    exp.txt <- FALSE
  }
  try.location <- try(match.arg(location, c("bottomleft", "topleft", 
                                            "bottomright", "topright")), silent = T)
  if (inherits(try.location, "try-error")) {
    stop(location, " is not an available named position to put the legend: please specify one of the following: 'bottomleft', 'topleft', 'bottomright', or 'topright'")
  }
  location <- try.location
  gen.args <- graphics:::.Pars
  poly.args <- names(formals(polygon))
  poly.args <- poly.args[-which(poly.args == "...")]
  text.args <- names(formals(text.default))
  text.args <- text.args[-which(text.args == "...")]
  if (is.null(breaks)) {
    if (is.null(val.range)) {
      val.range <- range(sim[[color.element]], na.rm = TRUE)
    }
    colramp <- colorRampPalette(col, alpha = T)(res)
    colramp <- alter.cols(colramp, alpha = evorates:::.lin.interp(alpha, 
                                                                  length(colramp)))
  }else {
    colramp <- colorRampPalette(col, alpha = T)(length(breaks) + 
                                                  1)
    colramp <- alter.cols(colramp, alpha = .lin.interp(alpha, 
                                                       length(colramp)))
    if (is.null(select.levels)) {
      select.levels <- 1:length(colramp)
    }
    select.levels <- select.levels[select.levels >= 1 & select.levels <= 
                                     (length(breaks) + 1)]
    select.levels <- sort(select.levels)
    colramp <- colramp[select.levels]
  }
  bds <- par("usr")
  bds.dims <- c(diff(bds[1:2]), diff(bds[3:4]))
  if (length(box.dims) == 0) {
    box.dims <- rep(NA, 2)
  }else if (length(box.dims) == 1) {
    box.dims <- c(box.dims, NA)
  }
  box.dims <- box.scale * ifelse(is.na(box.dims), bds.dims/c(30, 
                                                             5), box.dims)
  if (length(box.offset) == 0) {
    box.offset <- rep(NA, 2)
  }else if (length(box.offset) == 1) {
    box.offset <- c(box.offset, NA)
  }
  box.offset <- ifelse(is.na(box.offset), bds.dims/c(8, 20), 
                       box.offset)
  x.offset <- box.offset[1]
  y.offset <- box.offset[2]
  if (grepl("right", location)) {
    x.offset <- bds.dims[1] - box.dims[1] - box.offset[1]
  }
  if (grepl("top", location)) {
    y.offset <- bds.dims[2] - box.dims[2] - box.offset[2]
  }
  coords <- list(x = c(0, box.dims[1], box.dims[1], 0) + x.offset + 
                   bds[1], y = c(0, 0, box.dims[2], box.dims[2]) + y.offset + 
                   bds[3])
  y.int <- seq(coords$y[2], coords$y[3], length.out = length(colramp) + 
                 1)
  for (i in 1:length(colramp)) {
    do.call(polygon, c(x = list(coords$x), y = list(c(y.int[i], 
                                                      y.int[i], y.int[i + 1], y.int[i + 1])), border = NA, 
                       col = colramp[i], list(...)[!(names(list(...)) %in% 
                                                       c("x", "y", "border", "col", "adj")) & names(list(...)) %in% 
                                                     gen.args]))
  }
  do.call(polygon, c(x = list(coords$x), y = list(coords$y), 
                     col = NA, list(...)[!(names(list(...)) %in% c("x", "y", 
                                                                   "col", "adj")) & names(list(...)) %in% c(gen.args, 
                                                                                                            poly.args)]))
  side <- NA
  if (hasArg(side)) {
    if (list(...)$side <= 2 & list(...)$side >= 1) {
      side <- list(...)$side
    }
  }
  if (is.na(side)) {
    side <- if (grepl("left", location)) 
      2
    else 1
  }
  txt.args <- list(...)[!(names(list(...)) %in% c("x", "y", 
                                                  "labels")) & names(list(...)) %in% c(gen.args, text.args)]
  if (is.null(txt.args$adj) & is.null(txt.args$pos)) {
    txt.args$pos <- side * 2
  }
  if (!is.null(txt.col)) {
    txt.args$col <- txt.col
  }
  if (is.null(breaks)) {
    if (val.range[2] - val.range[1] == 0) {
      labels <- val.range[1]
    }else {
      if(is.na(labels_n)){
        labels <- pretty(seq(val.range[1], val.range[2], 
                             length.out = 100))
        labels <- labels[2:(length(labels) - 1)]
        
      }else{
        
        labels <- pretty(seq(val.range[1], val.range[2], 
                             length.out = 100), n=labels_n)
        labels <- labels[2:(length(labels) - 1)]
        
      }
    }
    y.pos <- coords$y[2] + (labels - val.range[1])/(diff(val.range)) * 
      (coords$y[3] - coords$y[2])
    if (exp.txt) {
      labels <- format(exp(labels), digits = 1)
    }
    do.call(text, c(x = coords$x[side], y = list(y.pos), 
                    labels = list(labels), txt.args))
  }else {
    labels <- paste(signif(breaks[-length(breaks)], 3), signif(breaks[-1], 
                                                               3), sep = " - ")
    labels <- c(paste("<", signif(breaks[1], 3)), labels, 
                paste(">", signif(breaks[length(breaks)], 3)))
    labels <- labels[select.levels]
    y.pos <- apply(cbind(y.int[-length(y.int)], y.int[-1]), 
                   1, mean)
    do.call(text, c(x = coords$x[side], y = list(y.pos), 
                    labels = list(labels), txt.args))
  }
  if (is.null(main)) {
    if (color.element == "R") {
      if (exp | exp.txt) {
        main <- expression(sigma^2)
      }else {
        main <- expression(ln ~ (sigma^2))
      }
    }else {
      main <- substitute(color.element)
    }
  }
  main.side <- NA
  if (hasArg(main.side)) {
    if (list(...)$main.side <= 4 & list(...)$main.side >= 
        1) {
      main.side <- list(...)$main.side
    }
  }
  if (is.na(main.side)) {
    main.side <- if (grepl("left", location)){ 
      2
    }else{ 4
    }
  }
  main.args <- list(...)[!(names(list(...)) %in% paste0("main.", 
                                                        c("x", "y", "labels", "side"))) & names(list(...)) %in% 
                           paste0("main.", c(gen.args, text.args))]
  names(main.args) <- gsub("main.", "", names(main.args))
  txt.args <- txt.args[!names(txt.args) %in% c("srt", "pos", 
                                               "adj", "offset")]
  txt.args[names(main.args) %in% names(txt.args)] <- main.args[names(main.args) %in% 
                                                                 names(txt.args)]
  main.args <- c(txt.args, main.args[!(names(main.args) %in% 
                                         names(txt.args))])
  if (is.null(main.args$adj) & is.null(main.args$pos)) {
    main.args$pos <- main.side
    if (main.side %in% c(2, 4) & is.null(main.args$srt)) {
      main.args$srt <- 90
      main.args$pos <- NULL
      main.args$adj <- c(0.5, c(-0.3, 1.3)[main.side/2])
    }
  }
  if (main.side %in% c(2, 4) & is.null(main.args$srt)) {
    main.args$srt <- 90
  }
  if (is.null(main.args$cex)) {
    main.args$cex <- 1
  }
  if (main.side %in% c(2, 4)) {
    x.pos <- coords$x[main.side/2]
    y.pos <- mean(coords$y[2:3])
    if (main.side == 4) {
      x.pos <- x.pos + bds.dims[1]/50
    }else {
      x.pos <- x.pos - bds.dims[1]/50
    }
  }else {
    x.pos <- mean(coords$x[1:2])
    y.pos <- coords$y[main.side/2 + 1.5]
    if (main.side == 3) {
      y.pos <- y.pos + bds.dims[2]/50
    }else {
      y.pos <- y.pos - bds.dims[2]/50
    }
  }
  do.call(text, as.list(c(x = x.pos, y = y.pos, labels = list(main), 
                          main.args)))
  invisible(coords)
}



plot_BePhyNE_simmap = function(contsimmaps, scale_atr=NA){
  #par(mfrow = c(2,2))
  
  char_names= colnames(contsimmaps[1,,1])
  char_space= list()
  count=0
  for (i in seq(1, length(char_names ), by = 2)) {
    count=count+1
    pred_chars = char_names[i:(i+1)]
    
    opt_BM= unlist(contsimmaps[,pred_chars[[1]],])
    brdth_BM = unlist(contsimmaps[,pred_chars[[2]],])
    
    traits_nichespace = do.call(rbind, lapply(1:length(opt_BM), function(i) backTransform1(cbind(opt_BM, brdth_BM, 1.0)[i,])))
    opt = traits_nichespace[,1]
    brdth = traits_nichespace[,2]
    
    if( sum(is.na(scale_atr))==1){
      char_space[[pred_chars[[1]] ]] = (opt  )
      char_space[[pred_chars[[2]] ]] = (brdth)
      
    }else{
      char_space[[pred_chars[[1]] ]] = (opt  * scale_atr$center[[count]]) + scale_atr$scale[[count]]
      char_space[[pred_chars[[2]] ]] = (brdth* scale_atr$scale[[count]])
    }
  }
  
  
  
  
  col_gradient = hcl.colors(100)
  for(i in 1:length(char_names)){
    
    contsimmap:::plot.contsimmap(traits = "phylogram", contsimmap = contsimmaps , Col.by = char_names[[i]], main = char_names[[i]], col=col_gradient,alpha=1)
    #contsimmap:::plot.contsimmap(contsimmap = sims, Col.by = char_names[[i]])
    
    legend.evorates_mod(sim=list(" " = char_space[[i]]), color.element = " ",col =col_gradient, location = "topleft",  exp.txt = F, labels_n=4)
    
    #col=<your color vector>,
    #breaks=<your custom vector of breaks, if you have one, otherwise it just picks 100 equally-spaced breaks like plot.contsimmap>)
    
  }
  
}
