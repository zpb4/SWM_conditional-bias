#function to create timeseries of different lags (- lag values) or leads (+ lag values) given an input timeseries
ts_lag_function<-function(ts,lag){
  ts_out<-c()
  if(lag<0){ts_out<-c(rep(0,abs(lag)),ts[1:(length(ts)+lag)])}
  if(lag>0){ts_out<-c(ts[(lag+1):length(ts)],rep(0,abs(lag)))}
  return(ts_out)
}

#function to create day-of-water-year (dowy) sequence given start and end dates
dowy_function<-function(start_date,end_date){
  idx<-seq(as.Date(start_date),as.Date(end_date),by='day')
  ixx<-as.POSIXlt(idx)
  
  ref_idx<-seq(as.Date('1983-10-01'),as.Date('1984-09-30'),by='day')
  ref_ixx<-as.POSIXlt(ref_idx)
  
  dowy_idx<-cbind(ref_ixx$mo,ref_ixx$mday,1:length(ref_idx))
  
  dowy_out<-rep(0,length(idx))
  
  for(i in 1:length(dowy_out)){
    inp_idx<-which(dowy_idx[,1]==ixx$mo[i] & dowy_idx[,2]==ixx$mday[i])
    dowy_out[i]<-dowy_idx[inp_idx,3]
  }
  return(dowy_out)
}

#AR simulation function
ar_sim<-function(n,ar_lags,ar_coef,innov){
  sim_vec<-rep(0,n+ar_lags)
  sim_vec[1:ar_lags]<-innov[sample(1:n,ar_lags)]
  for(i in 1:n){
    sim_vec[i+ar_lags]=t(sim_vec[(i+ar_lags-1):(i)])%*%ar_coef + innov[i]
  }
  return(sim_vec[(ar_lags+1):length(sim_vec)])
}