
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
