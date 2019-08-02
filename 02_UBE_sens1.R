# 02_UBE_sens1.R
# MR, May 2019

# Calibration and first step of the sensitivity analysis
# Script to be loaded on Ubelix to run simulation on different parameter sets

############################## For the UBELIX .job job
# Args from array in bash
args=(commandArgs(TRUE))
args=as.numeric(unlist(args))

print(args)

ps <- args[1] # parameter set

print(ps)

##############################
#directory to take the external data from (in ubelix it is in the same folder)
# dirl <- "tab_params/"
dirl <- ""
## requires paramters 
source("01_params_script_baselineparam.R")
#and the model script
source("03_functions_mainmodel.R")
# and the sensitivity parameter sets matrix
load(paste(dirl, "paramsets2019_01_26.RData", sep="") ) 

ptv <- c("alpha", "gamma_f", "gamma_m",  "omega", "tau", "pri1", "pri2", "pr23" , "pr2s", "pr3s", "prlcrc", "prrcdc", 
         "z1" ,"z2", "z3" ,"pi1", "pi2","pi3", "um1","um2","um3","dm1","dm2","dm3",
         "screening1", "screening2", "screening3", "screening4", "screening5","screening6", "screening7",
         "sensitivity", "specificity", "ve")

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

#loop over all HPV types

betaopt <- rep(NA, length(hpvtype))
pr3copt <- rep(NA, length(hpvtype))

for(ht in 1:length(hpvtype)){
  print(hpvtype[ht])
  ######################## Set the parameter right
  alpha <-  paramsets[ps,which(ptv=="alpha"),]
  gamma[1,] <- paramsets[ps,which(ptv=="gamma_f"),]
  gamma[2,] <- paramsets[ps,which(ptv=="gamma_m"),]
  omega[1,] <- paramsets[ps,which(ptv=="omega"),]
  omega[2,] <- paramsets[ps,which(ptv=="omega"),]
  tau[,1] <-  paramsets[ps,which(ptv=="tau"),]
  pri1[1,,1]  <-  paramsets[ps,which(ptv=="pri1"),]
  pri2[1,,1]  <-  paramsets[ps,which(ptv=="pri2"),]
  pr23[,1] <-   paramsets[ps,which(ptv=="pr23"),]
  pr2s[,1] <-  paramsets[ps,which(ptv=="pr2s"),]
  pr3s[,1] <-  paramsets[ps,which(ptv=="pr3s"),]
  prlcrc <-  paramsets[ps,which(ptv=="prlcrc"),1]
  prrcdc <-  paramsets[ps,which(ptv=="prrcdc"),1]
  z[1,1] <- paramsets[ps,which(ptv=="z1"),1]
  z[2,1] <- paramsets[ps,which(ptv=="z2"),1]
  z[3,1] <- paramsets[ps,which(ptv=="z3"),1]
  pi[1,1] <- paramsets[ps,which(ptv=="pi1"),1]
  pi[2,1] <- paramsets[ps,which(ptv=="pi2"),1]
  pi[3,1] <- paramsets[ps,which(ptv=="pi3"),1]
  um[1] <-  paramsets[ps,which(ptv=="um1"),1]
  um[2] <-  paramsets[ps,which(ptv=="um2"),1]
  um[3] <-  paramsets[ps,which(ptv=="um3"),1]
  dm[1] <-  paramsets[ps,which(ptv=="dm1"),1]
  dm[2] <-  paramsets[ps,which(ptv=="dm2"),1]
  dm[3] <-  paramsets[ps,which(ptv=="dm3"),1]
  s[1,1,1] <- paramsets[ps,which(ptv=="screening1"),1]
  s[1,2,1] <- paramsets[ps,which(ptv=="screening2"),1]
  s[1,3,1] <- paramsets[ps,which(ptv=="screening3"),1]
  s[1,4,1] <- paramsets[ps,which(ptv=="screening4"),1]
  s[1,5,1] <- paramsets[ps,which(ptv=="screening5"),1]
  s[1,6,1] <- paramsets[ps,which(ptv=="screening6"),1]
  s[1,7,1] <- paramsets[ps,which(ptv=="screening7"),1]
  
  ssens <-  paramsets[ps,which(ptv=="sensitivity"),1]
  sspec <-  paramsets[ps,which(ptv=="specificity"),1]
  ve    <-   paramsets[ps,which(ptv=="ve"),]
  
  
  #####################################
  vbeta <-  seq(0.5,1,0.1) #first step: says the range in which beta can be
  vpr3c <- c(0.0001, seq(0.007, 0.05,0.01),0.1, 0.25 ) #second step range of pr3c
  
  w <- ht
  y <-  which(dis_outcome=="cervix")
  prev_ipag2<- rep(NA, length(vbeta))
  inc_cc <- c()
  print(ht)
  inc<- c()
  
  ############################################ Calculates, for each individual HPV type at time, model outputs when varying either beta or pr3c
  # pr2c is considered to be 0.2*pr3c. Gives, the HPV prevalence for beta and cervical cancer incidence for pr3c, the linear regression between the two closest points estimates,
  #says if the rounded resulting estimates is the same as the data point estimate, indicates the lowest value possible (= 0), and the correct estimate for the data point to be fitted on.
  
  ########################### for varying beta
  # First check if with maximum beta (beta=1), the prevalence is higher than the one expected
  
  betap[1,w] <- vbeta[length(vbeta)]
  betap[2,w] <- vbeta[length(vbeta)]
  params <- c(w= ht, y = y )
  RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  ) 
  
  init <- RS$y
  spw_esp <- array(0, dim=c(  length(agegroups), length(comp)  ),
                   dimnames =list(  agegroups, comp) )
  
  for(i in 1: length(agegroups)) {
    for(j in 1: length(comp) ){
      spw_esp[i,j] <- ( (init[(i-1)*(length(comp)*4) + j] + init[(i-1)*(length(comp)*4) + (length(comp)*2)+j] ) / (sr*agr[i]) )* sagw[1+i]
    } }
  
  maxp <- spw_esp[3,2] / sagw[4]
  
  if(maxp < prev_pt_2024[ht] ){
    stop <- FAIL
  }else 
  {
    prev_ipag2[length(vbeta)] <- spw_esp[3,2] / sagw[4]
  }
  
  for(l in 1:(length(vbeta)-1)){
    betap[1,w] <- vbeta[l]
    betap[2,w] <- vbeta[l]
    params <- c(w= ht, y = y )
    RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  ) 
    init <- RS$y
    spw_esp <- array(0, dim=c(  length(agegroups), length(comp)  ),
                     dimnames =list(  agegroups, comp) )
    
    for(i in 1: length(agegroups)) {
      for(j in 1: length(comp) ){
        spw_esp[i,j] <- ( (init[(i-1)*(length(comp)*4) + j] + init[(i-1)*(length(comp)*4) + (length(comp)*2)+j] ) / (sr*agr[i]) )* sagw[1+i]
      } }
    
    prev_ipag2[l] <- spw_esp[3,2] / sagw[4]
  } 
  
  
  clos_smal <- which.min( (prev_pt_2024[ht] - prev_ipag2 )[if(prev_pt_2024[ht]< prev_ipag2[1]){1} else {( prev_pt_2024[ht] - prev_ipag2 )>0}])  #closest point to the given prev which is smaller than prev (or the closest point if there is no point smaller than prev)
  lr <- lm(prev_ipag2[c( clos_smal, clos_smal+1) ]~ vbeta[ c( clos_smal, clos_smal+1) ])  #linear regression between the two closest points estimates
  lincorr <- round(prev_ipag2[3],2) ==  round(lr$coefficients[2] * vbeta[3] + lr$coefficients[1],2)
  minbeta <- -lr$coefficients[1] /lr$coefficients[2]
  betadp <- (-lr$coefficients[1]+prev_pt_2024[ht] ) /lr$coefficients[2]
  prev <- prev_ipag2
  
  print(betadp)
  print(prev)
  
  betadp <- as.numeric(betadp)
  betadp
  
  ########################### for varying pr3c and pr2c  using scaled value for beta
  init_a <- array(0, dim=c(length(vpr3c),length(col_nam[-1]) ), dimnames = list(vpr3c,col_nam[-1]) )
  
  if(ht == 11 | ht == 12  )
  {
    pr3cdp <- 0
  }else{
    for(l in 1:length(vpr3c)){
      betap[1,w] <- betadp
      betap[2,w] <- betadp
      pr3c[ht] =vpr3c[l]
      pr2c[ht] <- pr3c[ht]*0.2
      params <- c(w= ht, y = y ) 
      RS <- runsteady(y = I0, fun = HPV_1,  parms = params, times = c(0, 1e5)  ) 
      init <- RS$y
      spw_esp <- array(0, dim=c(  length(agegroups), length(comp)  ),
                       dimnames =list(  agegroups, comp) )
      
      for(i in 1: length(agegroups)) {
        for(j in 1: length(comp) ){
          spw_esp[i,j] <- ( (init[(i-1)*(length(comp)*4) + j] + init[(i-1)*(length(comp)*4) + (length(comp)*2)+j] ) / (sr*agr[i]) )* sagw[1+i]
        } }
      inc_cc[l] <- ( ( z[1,y]*sum(spw_esp[,5]) + z[2,y]*sum(spw_esp[,7]) + z[3,y]*sum(spw_esp[,9] ) ) * 1e5 )
    }
    dp <- cc_inc_pertype_1099[ht]
    clos_smal <- which.min((dp - inc_cc )[(dp - inc_cc)>0])  #closest point to the given prev which is smaller than prev
    lr <- lm(inc_cc[c( clos_smal, clos_smal+1) ]~ vpr3c[ c( clos_smal, clos_smal+1) ]) 
    lincorr <- round(inc_cc[3],2) ==  round(lr$coefficients[2] * vpr3c[3] + lr$coefficients[1],2)
    minbeta <- -lr$coefficients[1] /lr$coefficients[2]
    pr3cdp <- (-lr$coefficients[1]+dp ) /lr$coefficients[2] 
    inc <- inc_cc
  }
  betaopt[ht] <- betadp
  pr3copt[ht] <- pr3cdp
}

tosave <- list(betaopt, pr3copt ) 
save( tosave , file="p.RData")




