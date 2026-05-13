

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
  
  X= lapply(1:nrow(traits[[1]]), function(sp)  lapply(3:length(pa_data[[sp]]), function(pred) pa_data[[sp]][[pred]]))
  
  betas<-lapply(1:length(traits), function(pred) traits2coefs(traits[[pred]]))
  
  ###glm stuff
  yy_prelogit_sep  = lapply( 1:nrow(traits[[1]]), function(sp)  lapply(1:length(traits), function(pred) (betas[[pred]][sp,1] + betas[[pred]][sp,2]*X[[sp]][[pred]]+ betas[[pred]][sp,3]*(X[[sp]][[pred]]^2)  ) ) )
  
  yy_prelogit= lapply(1:nrow(traits[[1]]), function(sp) rowSums(matrix(unlist(yy_prelogit_sep[[sp]]), ncol=length(yy_prelogit_sep[[sp]]), byrow=F)))
  
  presProb <-lapply(1:nrow(traits[[1]]), function(sp) data.frame(y=pa_data[[sp]]$y, fitted.values=1/(1+exp(-1*yy_prelogit[[sp]] ) )) ) #convert to probability using logit link
  ###
  
  return(presProb)
  
}


assess.predict = function(presProb, plot=F){
  
  if(sum(is.na(presProb$y))==1){
    
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


predict_stats  = function(traits, pa_data){
  
  sp_names = unlist(lapply(pa_data, function(i) i$species))
  
  presProb = predict.ENE(traits, pa_data)
  
  predict_list =  lapply(1:length(presProb), function(i) {
    #print(i)
    assess.predict(presProb[[i]])
  }
  )
  
  names(predict_list) = sp_names
  
  return(predict_list)
  
}


AUC_posterior_median =function(log_summary, pa_data){
  
  traits = lapply(log_summary$median_parlist$traits, function(i) i[[1]])
  
  predict_stats_list = predict_stats(traits, pa_data)
  
  return(predict_stats_list)
}



