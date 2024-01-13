#load required packages
library(stringr)
library(ranger)

#set working directory
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
#q_obs[q_obs==-9999]<-0
q_obs[q_obs<=0]<-min(q_obs[q_obs>0])
q_sim[q_sim<=0]<-min(q_sim[q_sim>0])

#---------------------------------------------------------
#1) Fit conditional bias model to historical Qobs and Qsim
cbias_out<-fit_cbias(Qobs = q_obs,Qsim = q_sim,start_date = start_date,end_date = end_date)

cbias_fit<-cbias_out[[1]] #fitted model to use for 'sim_cbias' function demonstrated below
cbias_qsim<-cbias_out[[2]] #estimated conditionally debiased simulation value

#---------------------------------------------------------
#2) Simulate from scenario folder
#code below works if you extract the entire scenario folder from the Box repository ('SAC-SMA_WGENext_15CDEC')
scenario<-1 #select the appropriate scenario
Q_sim_scenario<-read.table(paste('data/',scenario,'/simflow_sacsma_SHA.txt',sep='')) #read in Qsim data

#extract dates from Qsim_scenario data:
qsim_scenario_st_date<-paste(Q_sim_scenario[1,1],'-',str_pad(Q_sim_scenario[1,2],2,pad='0'),'-',str_pad(Q_sim_scenario[1,3],2,pad='0'),sep='')
qsim_scenario_end_date<-paste(Q_sim_scenario[dim(Q_sim_scenario)[1],1],'-',str_pad(Q_sim_scenario[dim(Q_sim_scenario)[1],2],2,pad='0'),'-',str_pad(Q_sim_scenario[dim(Q_sim_scenario)[1],3],2,pad='0'),sep='')

#simulate the conditionally debiased Qsim with the 'sim_cbias' function
st=qsim_scenario_st_date
end=qsim_scenario_end_date

#extract flows from column 4
Q_sim_scenario_flow<-Q_sim_scenario[,4]
#correct any zero flows or negative values to minimum non-zero flow
Q_sim_scenario_flow[Q_sim_scenario_flow<=0]<-min(Q_sim_scenario_flow[Q_sim_scenario_flow>0])

#simulate new conditionally debiased data
Q_sim_scenario_cbias<-sim_cbias(cbias_fit = cbias_fit, Qsim = Q_sim_scenario_flow, start_date = st, end_date = end)

#----------------------------------------------------------------------------------------
#3) Save to a .csv file for compatibility with other programs
#output .csv file of conditionally debiased simulations
dates<-as.character(seq(as.Date(st),as.Date(end),by='day'))
data_out<-cbind(dates,Q_sim_scenario_flow,Q_sim_scenario_cbias)
colnames(data_out)<-c('date','Qsim','Qsim_cbias')

#save to 'out' repo with scenario name appended to .csv title
write.csv(data_out,paste('./out/SWM-data_cond-debias_scenario-',scenario,'.csv',sep=''),row.names = FALSE)


#--------------------------------------------------END------------------------------------------