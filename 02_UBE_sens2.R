# 02_UBE_sens2.R
# MR, May 2019

# Script to be loaded on Ubelix to run simulation on different parameter sets
# Uses the parameters sets which achieve calibration to load the vacc scenarios and costs sensitivity analysis

############################## For the UBELIX .job job
# Args from array in bash
args=(commandArgs(TRUE))
args=as.numeric(unlist(args))

print(args)
ps <- args[1] # parameter set (from working parameterset)

print(ps)
##############################
#directory to take the external data from (in ubelix it is in the same folder)
# dirl <- "tab_params/"
dirl <- ""
## requires paramters 
source("01_params_script_baselineparam.R")
#and the model script
source("03_functions_mainmodel.R")
#function for costs (to groups the different hpvtypes)
#source("03_functions_costs.R") #load it later

load(paste(dirl, "paramsets_ok2019_06_04.RData", sep="") )
load(paste(dirl, "beta_ok2019_06_04.RData", sep="") )
load(paste(dirl, "pr3c_ok2019_06_04.RData", sep="") ) 

# dim(paramsets_ok)
print(ps)

dim(beta_ok)
dim(pr3c_ok)

ptv <- c("alpha", "gamma_f", "gamma_m",  "omega", "tau", "pri1", "pri2", "pr23" , "pr2s", "pr3s", "prlcrc", "prrcdc", 
         "z1" ,"z2", "z3" ,"pi1", "pi2","pi3", "um1","um2","um3","dm1","dm2","dm3",
         "screening1", "screening2", "screening3", "screening4", "screening5","screening6", 
         "screening7", #withdraw this if not vtest
         "sensitivity", "specificity", "ve")

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

#define different vaccination scenarios:
diffv <- c(0.56, 0.8)
vaccstrat <-  c( paste( rep("V4v", each =length(diffv) ), diffv*100, sep="_"), paste( rep("V9v", each =length(diffv) ), diffv*100, sep="_") )

#Set time after vaccine introduction
time <- seq(0, 200, 1)

#function for costs (to groups the different hpvtypes)
source("03_functions_costs.R")

#store data
tsave1 <- array(0, dim=c(length(time),length(col_nam),length(vaccstrat), length(hpvtype)), dimnames = list(time,col_nam, vaccstrat, hpvtype) )  

######################## Set the parameter right
# ptv
alpha <-  paramsets_ok[ps,which(ptv=="alpha"),]
gamma[1,] <- paramsets_ok[ps,which(ptv=="gamma_f"),]
gamma[2,] <- paramsets_ok[ps,which(ptv=="gamma_m"),]
omega[1,] <- paramsets_ok[ps,which(ptv=="omega"),]
omega[2,] <- paramsets_ok[ps,which(ptv=="omega"),]
tau[,1] <-  paramsets_ok[ps,which(ptv=="tau"),]
pri1[1,,1]  <-  paramsets_ok[ps,which(ptv=="pri1"),]
pri2[1,,1]  <-  paramsets_ok[ps,which(ptv=="pri2"),]
pr23[,1] <-   paramsets_ok[ps,which(ptv=="pr23"),]
pr2s[,1] <-  paramsets_ok[ps,which(ptv=="pr2s"),]
pr3s[,1] <-  paramsets_ok[ps,which(ptv=="pr3s"),]
prlcrc <-  paramsets_ok[ps,which(ptv=="prlcrc"),1]
prrcdc <-  paramsets_ok[ps,which(ptv=="prrcdc"),1]
z[1,1] <- paramsets_ok[ps,which(ptv=="z1"),1]
z[2,1] <- paramsets_ok[ps,which(ptv=="z2"),1]
z[3,1] <- paramsets_ok[ps,which(ptv=="z3"),1]
pi[1,1] <- paramsets_ok[ps,which(ptv=="pi1"),1]
pi[2,1] <- paramsets_ok[ps,which(ptv=="pi2"),1]
pi[3,1] <- paramsets_ok[ps,which(ptv=="pi3"),1]
um[1] <-  paramsets_ok[ps,which(ptv=="um1"),1]
um[2] <-  paramsets_ok[ps,which(ptv=="um2"),1]
um[3] <-  paramsets_ok[ps,which(ptv=="um3"),1]
dm[1] <-  paramsets_ok[ps,which(ptv=="dm1"),1]
dm[2] <-  paramsets_ok[ps,which(ptv=="dm2"),1]
dm[3] <-  paramsets_ok[ps,which(ptv=="dm3"),1]
s[1,1,1] <- paramsets_ok[ps,which(ptv=="screening1"),1]
s[1,2,1] <- paramsets_ok[ps,which(ptv=="screening2"),1]
s[1,3,1] <- paramsets_ok[ps,which(ptv=="screening3"),1]
s[1,4,1] <- paramsets_ok[ps,which(ptv=="screening4"),1]
s[1,5,1] <- paramsets_ok[ps,which(ptv=="screening5"),1]
s[1,6,1] <- paramsets_ok[ps,which(ptv=="screening6"),1]
s[1,7,1] <- paramsets_ok[ps,which(ptv=="screening7"),1]#withdraw this if not vtest

ssens <-  paramsets_ok[ps,which(ptv=="sensitivity"),1]
sspec <-  paramsets_ok[ps,which(ptv=="specificity"),1]
ve    <-   paramsets_ok[ps,which(ptv=="ve"),]

betap[1,] <- beta_ok[,ps]
betap[2,] <- beta_ok[,ps]

pr3c <- pr3c_ok[,ps]
pr2c <- pr3c*0.2

pb <- txtProgressBar(min = 0, max = length(hpvtype), style = 3)

#loop over all HPV types
for(htt in 1:(length(hpvtype)-2)){
  ht <- c(1:13)[-c(11,12)][htt]
  print(hpvtype[ht])
  ######################################### STEADY state for all HPV types, before starting preventions intervention (vaccination)
  #set vaccination to 0 (here for quadrivalent vaccination, all other already are at 0)
  #start with setting vaccination to 0 
  pv <- array(0, dim=c(length(sex),length(hpvtype), length(vaccines)), dimnames = list(sex,hpvtype, vaccines )  )
  vs <- 1
  w <- ht
  y <- 1
  
  params <- c(w= w, y = y )
  #this runs it until steady state and gives the initial vector at steady state takes about 2 min
  RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  )   
  #set the new initial values based on steady-state
  init <- RS$y
  
  tempin <- c()
  for(k in 0:(length(agegroups)-1) ){
    tempin[k+1] <- sum(init[ c(k*(13*4)+5,k*(13*4)+7 ,k*(13*4)+9, k*(13*4)+(13*2)+5,
                               k*(13*4)+(13*2)+7 ,k*(13*4)+(13*2)+9)]   * rep( z[,y], 2) / (sr*agr[1+k])  * sagw[2+k])
  }
  
  print(sum(tempin))
  
  # dperHPVt <- array(0, dim=c(length(time),length(col_nam),length(vaccstrat)), dimnames = list( time,col_nam,  vaccstrat)  )
  # First vaccination scenarios with V4v vaccination
  for (j in 1:length(diffv) ){
    w <- ht
    y <- 1
    params <- c(w= w, y = y )
    vs <- 2
    pv[1,,vs] <- 0
    pv[1,c(1,2,11,12),vs] <- diffv[j]
    sim <- as.data.frame(ode(init, time, HPV_1, parms=params))
    # dperHPVt[,,j] <- as.matrix(sim)
    tsave1[,,j,ht] <- as.matrix(sim)
    #########################################nonavalent vacc
    vs <- 3
    pv[1,,vs] <- diffv[j]
    pv[1,c(5,6,8,13),vs] <- 0
    sim <- as.data.frame(ode(init, time, HPV_1, parms=params))
    # dperHPVt[,,length(diffv)+j] <- as.matrix(sim)
    tsave1[,,length(diffv)+j,ht] <- as.matrix(sim)
    
  }
  
  setTxtProgressBar(pb, ht)
  
}  

# dperHPVt[1,1:10,]
dim(tsave1)

tsave1[1,1:5,1,]

############### Vérifier et faire en sorte que on ait les compartiments dont on a besoin par la suite
#vector of element which are not SIR nor deaths, over sex, riskgroups and age groups



tsave2 <- array(NA, dim=c(length(time),length(col_nam2),length(vaccstrat)), dimnames = list(time,col_nam2, vaccstrat) )  

is.na(tsave2)
dim(tsave2)

tsave2 <- dfin_group_f(tsave1)  

dim(tsave2)
tsave2[2,1:10,1]
############################################

tosave <- list(tsave2 ) 

save( tosave , file="p.RData")

#load("dfinal2.Rdata")
# dim(temp2)  
# 
# temp1 <- dfinal
# temp2 <- dfinal2
# save( dfinal3 , file="dfinal3_12.RData")
# save( dfinal2 , file="dfinal2.RData")
# save( dfinal , file="dfinal.RData")
