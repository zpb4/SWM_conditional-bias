
source('./src/hybrid-loess-fit_functions.R')


#Function to estimate conditional expectation of observations given simulations

fit_cbias<-function(Qobs,Qsim,start_date,end_date){
  
  idx<-seq(as.Date(start_date),as.Date(end_date),by='day')
  ixx<-as.POSIXlt(idx)

  cbias_fit<-vector('list',12)
  cbias_vec<-Qsim

  #fit hybrid LOESS model to the data (Qsim v Qobs)
  for(i in 1:12){
    seas<-which(ixx$mon==(i-1))

    #try symmetric (non-parametric) LOESS fit
    hy_fit<-try(hyb_loess_fit(Qsim[seas],Qobs[seas],.75,2,'symmetric'),T)
    #if that doesn't work, try a gaussian fit
    if(class(hy_fit)=='try-error'){
      print(paste('mth',i,'trying gaussian'))
      hy_fit<-try(hyb_loess_fit(Qsim[seas],Qobs[seas],.75,2,'gaussian'),T)
    }
    #if that doesn't work reduce to a local linear model (degree 1) fit
    if(class(hy_fit)=='try-error'){
      print(paste('mth',i,'trying deg1'))
      hy_fit<-try(hyb_loess_fit(Qsim[seas],Qobs[seas],.75,1,'gaussian'),T)
    }
    cbias_fit[[i]]<-hy_fit
    cbias<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],Qsim[seas])
    cbias_vec[seas]<-cbias
  }

#set any NAs or negative condition expectation values to the observed value
na_idx<-which(is.na(cbias_vec)==T)
neg_idx<-which(cbias_vec<0)

cbias_vec[na_idx]<-Qsim[na_idx]
cbias_vec[neg_idx]<-Qsim[neg_idx]

return(list(cbias_fit,cbias_vec))
}

sim_cbias<-function(cbias_fit,Qsim,start_date,end_date){
  idx<-seq(as.Date(start_date),as.Date(end_date),by='day')
  ixx<-as.POSIXlt(idx)
  
  cbias_vec<-Qsim
  
  #fit hybrid LOESS model to the data (Qsim v Qobs)
  for(i in 1:12){
    seas<-which(ixx$mon==(i-1))

    hy_fit<-cbias_fit[[i]]
    cbias<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],Qsim[seas])
    cbias_vec[seas]<-cbias
  }
  
  #set any NAs or negative condition expectation values to the observed value
  na_idx<-which(is.na(cbias_vec)==T)
  neg_idx<-which(cbias_vec<0)
  
  cbias_vec[na_idx]<-Qsim[na_idx]
  cbias_vec[neg_idx]<-Qsim[neg_idx]
  
  return(cbias_vec)
}


###########################################END################################