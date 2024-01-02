#Fit VAR and SGED models to Differenced errors
library(stringr)
library(ranger)

setwd('d:/SWM_conditional-bias/')

#read in data
Q_obs<-read.table('./data/FNF_SHA_mm.txt')
Q_sim<-read.table('./data/simflow_sacsma_SHA_short.txt')

#required functions
source('./src/cbias-hsked_functions.R')
source('./src/helper_functions.R')

#extract date indices from data files
qobs_st_date<-paste(Q_obs[1,1],'-',str_pad(Q_obs[1,2],2,pad='0'),'-',str_pad(Q_obs[1,3],2,pad='0'),sep='')
qobs_end_date<-paste(Q_obs[dim(Q_obs)[1],1],'-',str_pad(Q_obs[dim(Q_obs)[1],2],2,pad='0'),'-',str_pad(Q_obs[dim(Q_obs)[1],3],2,pad='0'),sep='')
qsim_st_date<-paste(Q_sim[1,1],'-',str_pad(Q_sim[1,2],2,pad='0'),'-',str_pad(Q_sim[1,3],2,pad='0'),sep='')
qsim_end_date<-paste(Q_sim[dim(Q_sim)[1],1],'-',str_pad(Q_sim[dim(Q_sim)[1],2],2,pad='0'),'-',str_pad(Q_sim[dim(Q_sim)[1],3],2,pad='0'),sep='')

qobs_idx<-seq(as.Date(qobs_st_date),as.Date(qobs_end_date),'day')
qsim_idx<-seq(as.Date(qsim_st_date),as.Date(qsim_end_date),'day')

#specify date range that is available in both Qobs and Qsim
start_date<-'1987-10-01'
end_date<-'2018-09-30'

#extract q_obs and q_sim vector for specified date range
q_obs<-Q_obs[which(qobs_idx==start_date):which(qobs_idx==end_date),4]
q_sim<-Q_sim[which(qsim_idx==start_date):which(qsim_idx==end_date),4]

#set any 0 or negative flow values to minimum non-zero flow
q_obs[q_obs<=0]<-min(q_obs[q_obs>0])
q_sim[q_sim<=0]<-min(q_sim[q_sim>0])

#Fit conditional bias model
cbias_out<-fit_cbias(Qobs = q_obs,Qsim = q_sim,start_date = start_date,end_date = end_date)

cbias_fit<-cbias_out[[1]] #fitted model 
cbias_qsim<-cbias_out[[2]] #estimated conditionally debiased simulation value

#calculate lambdas for both raw Qsim and debiased Qsim
lambda_raw<-log(q_sim/q_obs)
lambda_cbias<-log(cbias_qsim/q_obs)

#-------------------------------------------------------------------------------
#Random forest debiasing approach

#day-of-water-year (dowy) predictor variable
dowy_out<-dowy_function(start_date = start_date, end_date = end_date)

#Note: Tried RF prediction using dowy, lag -1:-3 and +1:+3 predictors; found that Qsim, lag -1, lag +1 contained most of the variable importance
#rf_predictor_mat<-cbind(q_sim,dowy_out,ts_lag_function(q_sim,-3),ts_lag_function(q_sim,-2),ts_lag_function(q_sim,-1),ts_lag_function(q_sim,1),ts_lag_function(q_sim,2),ts_lag_function(q_sim,3))
#colnames(rf_predictor_mat)<-c('Qsim','dowy','lag3','lag2','lag1','fwd1','fwd2','fwd3')

#reduced set of predictors
rf_predictor_mat<-cbind(q_sim,ts_lag_function(q_sim,-1),ts_lag_function(q_sim,1))
colnames(rf_predictor_mat)<-c('Qsim','lag1','fwd1')

#fit the Random Forest model between predictors and Qobs
rf_fit<-ranger(x=rf_predictor_mat,y=q_obs,importance = 'impurity')
#RF model predictions of Qobs (ie debiased Qsim)
rf_cbias_qsim<-predict(rf_fit,data=rf_predictor_mat)$predictions
#calculate lambdas with RF debiased simulations
lambda_rf_cbias<-log(rf_cbias_qsim/q_obs)

#-----------------------------------------------------------------------
#Plots

#log(q_obs) vs lambda plots
par(mfrow=c(3,1))
plot(log(q_obs),lambda_raw,ylim=c(-5,5))
plot(log(q_obs),lambda_cbias,ylim=c(-5,5))
plot(log(q_obs),lambda_rf_cbias,ylim=c(-5,5))

#lambda histogram plots
par(mfrow=c(3,1))
brks<-c(-10,seq(-5,5,0.25),10)

hist(lambda_raw,breaks=brks,xlim=c(-5,5))
hist(lambda_cbias,breaks=brks,xlim=c(-5,5))
hist(lambda_rf_cbias,breaks=brks,xlim=c(-5,5))

#plot conditionally debiased results comparison on smaller x-axis
par(mfrow=c(2,1))
brks<-c(-10,seq(-5,5,0.1),10)

hist(lambda_cbias,breaks=brks,xlim=c(-2,2))
hist(lambda_rf_cbias,breaks=brks,xlim=c(-2,2))


#plot variable importance for RF model
var_vec<-colnames(rf_predictor_mat)
var_imp<-rf_fit$variable.importance / sum(rf_fit$variable.importance)
srt_var_imp<-sort(var_imp,index.return=T,decreasing=T)

par(mfrow=c(1,1))
barplot(srt_var_imp$x,ylim=c(0,0.35),names.arg = var_vec[srt_var_imp$ix],main='',
        xlab='',ylab='Variable Importance (fraction)')


#----------------------------------------------------------------------------
#some testing of output

#fit AR(3) model to each lambda
arfit_raw<-arima(lambda_raw,order=c(3,0,0),include.mean = F)
arfit_cbias<-arima(lambda_cbias,order=c(3,0,0),include.mean = F)
arfit_rf_cbias<-arima(lambda_rf_cbias,order=c(3,0,0),include.mean = F)

#bootstrap resample residuals
raw_resids_resamp<-arfit_raw$residuals[sample(1:length(q_sim),length(q_sim),replace = TRUE)]
cbias_resids_resamp<-arfit_cbias$residuals[sample(1:length(q_sim),length(q_sim),replace = TRUE)]
rf_cbias_resids_resamp<-arfit_rf_cbias$residuals[sample(1:length(q_sim),length(q_sim),replace = TRUE)]

#simulate new lambda timeseries
lambda_raw_sim<-ar_sim(length(q_sim),3,arfit_raw$coef,raw_resids_resamp)
lambda_cbias_sim<-ar_sim(length(q_sim),3,arfit_cbias$coef,cbias_resids_resamp)
lambda_rf_cbias_sim<-ar_sim(length(q_sim),3,arfit_rf_cbias$coef,rf_cbias_resids_resamp)

#simulate from SWM
q_swm_raw = q_sim/exp(lambda_raw_sim)
q_swm_cbias = cbias_qsim/exp(lambda_cbias_sim)
q_swm_rf_cbias = rf_cbias_qsim/exp(lambda_rf_cbias_sim)

#compare maximum values
max(q_obs)
max(q_swm_raw) #very large, can be > 10X q_obs
max(q_swm_cbias)
max(q_swm_rf_cbias)

#mean values
mean(q_obs)
mean(q_swm_raw) #well above mean of q_obs, driven high by high maximum values
mean(q_swm_cbias) #pretty close
mean(q_swm_rf_cbias) #pretty close but low

#median values
median(q_obs)
median(q_swm_raw) #median is not bad
median(q_swm_cbias) #pretty close
median(q_swm_rf_cbias) #pretty close but lower


##################################################END################################################