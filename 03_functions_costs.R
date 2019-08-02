# 03_functions_costs.R
# MR, May 2019

########### For loading and showing costs: function adding HPV types and calculating costs


################## Pre functions and codes to load the function which gives the costs and qalys -----------------------------------------------
##################### 1st step is to combine the different compartements from dfinal across the different hpv types
##################### The disease stages from the different types will be added up. Vaccination uptake has to be accounted once, and the rest of the
##################### compartements, except the deaths, will be groupes (SIR) and calculated as 1 - (disease stages + vaccination uptake)
#function which gives the emplacement of specific compartements over age groupes sex risk groups. 
fagmf <- function(comptla, compartements){
  nbcompt <- rep(NA, length(agegroups)*length(riskgroups))
  for(i in 0: (length(agegroups)-1) ){
    nbcompt[ ( i*(length(comptla)*4) +1 ) : ( ( i*(length(comptla)*4) +1 )+ 3 ) ] <- 1+ c( (i*length(compartements)*4 + comptla), #female low
                                                                                           (i*length(compartements)*4 + length(compartements)*1 + comptla), #male low
                                                                                           (i*length(compartements)*4 + length(compartements)*2 + comptla), #female high
                                                                                           (i*length(compartements)*4 + length(compartements)*3 + comptla))  #male high
  }
  return(nbcompt)
}
# fagmf(1, length(newcomp) )

#same only for women
fag <- function(comptla, compartements){
  nbcompt <- rep(NA, length(agegroups)*length(riskgroups))
  for(i in 0: (length(agegroups)-1) ){
    nbcompt[ ( i*(length(comptla)*2) +1 ) : ( ( i*(length(comptla)*2) +1 )+ 1 ) ] <- 1+ c( (i*length(compartements)*4 + comptla),                                                                                        (i*length(compartements)*4 + length(compartements)*2 + comptla)) 
  }
  return(nbcompt)
}

# fag(3, newcomp)
# col_nam2[ c( fag(3, newcomp)) ]

########### new compartements: Groups SIR together, and disease stages remain and will be added for each types, (vaccination stays), deaths are excluded
newcomp <- c("SIR", comp[3:10], comp[13])
col_nam2 <- compnamef(newcomp) 

#vector of element which are not SIR nor deaths, over sex, riskgroups and age groups
vnsir <- array(0, dim=c(length(3:11), length(agegroups)*4 ) , dimnames = list(col_nam2[3:11], col_nam2[fagmf(2,newcomp)] ) ) #without the DCC (deaths from cc)
for(i in 1:9){
  vnsir[i,] <- fagmf(i+1, newcomp)
}
vnsir

## These for loops, take together, for each sex, riskgroups, agegroups, the disease stages for all different types, the % vaccinated and the number of deaths and makes 1- the sum
#of all this to have the SIR grouped, for the costs calculations. It groups it according to the newcomp vector

dfinal2 <- array(0, dim=c(length(time),length(col_nam2),length(vaccstrat)), dimnames = list(time,col_nam2, vaccstrat) )  

dfin_group_f <- function(dfinal){
  for(t in 1:length(time) ){
    for(i in 1:length(vaccstrat) ) {
      for( j in c(3:10 )){
        dfinal2[t,fagmf(j-1, newcomp),i] <-  rowSums(dfinal[t,fagmf(j, comp),i,] )  
      }# "CIN2" "CIN3" "ULCC" "DLCC" "URCC" "DRCC" "UDCC" "DDCC" "DCC" (comp[3:11])
      dfinal2[t,fagmf(13-3, newcomp),i] <- (dfinal[t,fagmf(13, comp),i,1] )  
      #vaccinated    
      for(l in 1:length(vnsir[1,]) ){
        dfinal2[t,fagmf(1, newcomp)[l],i] <-   as.vector(N)[l] - sum( dfinal2[t, as.vector(vnsir[,l]) ,i  ] ) #for sex, risk group and age groups, SIR is
      } } }  #for sex, risk group and age groups, SIR is N[sex, riskgroup, agegroup] - sum(CIN2, CIN3, ULCC, etc.)
  return(dfinal2) } #dfinal is an array of dim [time, length(col_nam), length(vaccstrat), length(hpvtype)]

# dfinal2[1,,1]
# sum(dfinal2[1,,1][fag(1, newcomp)] ) * popmod *1e6
# sum(dfinal2[200,,1][fag(11, newcomp)] ) * popmod *1e6
# dfinal <- tsave1

# dfin_group_f()

################### Function for load and show costs

#Big function giving, for women, cervical cancer and women vaccination the costs and icer
diffcosts <- c("pap_costs", "CIN2CIN3_costs", "cancer_costs", "vacc_costs", "AE_costs", "AE_QALY", "TOT_QALY", "TOT_costs")

cqoutcome <- array(0, dim=c(length(time),length(diffcosts), length(vaccstrat) ), dimnames = list(time,diffcosts, vaccstrat) )
icer  <- matrix(0, nrow= length(time), ncol= length(vaccstrat) )
icer2 <- matrix(0, nrow= length(time), ncol= length(vaccstrat) )
icer100 <- matrix(NA, ncol=length(vaccstrat), nrow=2); colnames(icer100) <- vaccstrat;rownames(icer100) <- c("50%V4v_baseline", "novacc_baseline")

#function that gives the costs for the different vaccination strategies or parameters, and the ICER. Only women vaccination. FOr men vaccination, another function is needed
# yp <- 1


################## Nouvelle version avec age weight ajusté

#### NF with more parameters to be defined
fcostsvaccstrat_varyp2 <- function(dfinal2, nbsens, C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q){
  for(i in 1:length(vaccstrat) ){
    pv <- array(0, dim=c(length(sex),length(hpvtype), length(vaccines)), dimnames = list(sex,hpvtype, vaccines )  )
    if(i == 1 | i == 2){
      vs <- 2
      pv[1,1,vs] <- diffv[i]
    }else if(i == 3| i == 4){
      vs <- 3
      pv[1,1,vs] <- diffv[i-2]
    }
    dataf <- dfinal2[,,i]
    # pap_costs
    cqoutcome[,1,i] <- apply(dataf, 1,
                             function(x) 
                               (sum( rep( rep(c(paramsets_ok[nbsens,which(ptv=="screening1"),1],
                                                paramsets_ok[nbsens,which(ptv=="screening2"),1],
                                                paramsets_ok[nbsens,which(ptv=="screening3"),1],
                                                paramsets_ok[nbsens,which(ptv=="screening4"),1],
                                                paramsets_ok[nbsens,which(ptv=="screening5"),1],
                                                paramsets_ok[nbsens,which(ptv=="screening6"),1],
                                                paramsets_ok[nbsens,which(ptv=="screening7"),1]
                               )/ agr * sagw[2:8], each=length(riskgroups)),2) 
                               *  Cpt*(1- paramsets_ok[nbsens,which(ptv=="sensitivity"),1])
                               * x[ c( fag(2, newcomp), fag(3, newcomp) ) ]  )   #costs of false negative screens  * ()
                               +  sum( rep( rep(c(paramsets_ok[nbsens,which(ptv=="screening1"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening2"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening3"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening4"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening5"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening6"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening7"),1]
                               ) / agr * sagw[2:8], each=length(riskgroups)),2) *  
                                 ( x[c( fag(1, newcomp), fag(10, newcomp) )] - 
                                     ( (1-paramsets_ok[nbsens,which(ptv=="specificity"),1]) 
                                       * x[c( fag(1, newcomp), fag(10, newcomp) )] ) )  *Cpt )  #for SIR and vaccinated, cost negative screens
                               +  sum( rep( rep(c(paramsets_ok[nbsens,which(ptv=="screening1"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening2"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening3"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening4"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening5"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening6"),1],
                                                  paramsets_ok[nbsens,which(ptv=="screening7"),1]
                               ) / agr * sagw[2:8], each=length(riskgroups)),2) *
                                 ( (1-paramsets_ok[nbsens,which(ptv=="sensitivity"),1]) 
                                   * x[c( fag(1, newcomp), fag(10, newcomp) )] )  * Cfp ) ) ) #for SIR and vaccinated, cost of follow up of false positives
    
    # CIN2CIN3_costs 
    cqoutcome[,2,i] <- apply(dataf, 1, function(x) (  sum( rep( rep(c(paramsets_ok[nbsens,which(ptv=="screening1"),1],
                                                                      paramsets_ok[nbsens,which(ptv=="screening2"),1],
                                                                      paramsets_ok[nbsens,which(ptv=="screening3"),1],
                                                                      paramsets_ok[nbsens,which(ptv=="screening4"),1],
                                                                      paramsets_ok[nbsens,which(ptv=="screening5"),1],
                                                                      paramsets_ok[nbsens,which(ptv=="screening6"),1],
                                                                      paramsets_ok[nbsens,which(ptv=="screening7"),1]
    )/ agr *sagw[2:8], each=2),2)*paramsets_ok[nbsens,which(ptv=="sensitivity"),1]* 
      rep( Ct[1:2], each=length(riskgroups)*length(agegroups) )* x[c( fag(2, newcomp), fag(3, newcomp))]))) #costs true positive treatement costs
    # cancer_costs 
    cqoutcome[,3,i] <- apply(dataf, 1, function(x) (  sum( rep( Ct[3:5,y], each=length(riskgroups)*length(agegroups) )*
                                                             rep( c(paramsets_ok[nbsens,which(ptv=="z1"),1],
                                                                    paramsets_ok[nbsens,which(ptv=="z2"),1],
                                                                    paramsets_ok[nbsens,which(ptv=="z3"),1]), each=length(riskgroups)*length(agegroups) )*
                                                             x[c( fag(4, newcomp), fag(6, newcomp), fag(8, newcomp)  )] * rep(rep(sagw[2:8] / agr, each=2 ),3) )) ) 
    # vacc_costs 
    cqoutcome[,4,i] <- apply(dataf, 1, function(x) ifelse(vs==2, sum(  C4v* (as.vector(pv[,1,vs]) *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp  )),#yp is an additional factor to adapt the number of school aged population to be vaccinated (reversed age-pyramid)
                                                          sum( C9v * (as.vector(pv[,1,vs]) *g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp  )) ) )
    # AE_costs 
    cqoutcome[,5,i]  <- apply(dataf, 1, function(x) sum( (rme*Cme + rse*Cse)* (as.vector(pv[,1,vs])*g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp )) ) 
    # AE_QALY 
    cqoutcome[,6,i] <- apply(dataf, 1, function(x)  sum( (rme*Qme + rse*Qse)* (as.vector(pv[,1,vs])*g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp ) )) 
    
    # TOT_QALY
    cqoutcome[,7,i] <- apply(dataf, 1, function(x) (  sum(x[c(fagmf(1, newcomp) ,fagmf(10, newcomp) ) ] * rep(rep(sagw[2:8]/ agr, each=4 ),2) ) #sum healthy women and men
                                                      + sum( rep( q[1:2,y], each= length(riskgroups)*length(agegroups) )* x[c(fag(2, newcomp), fag(3, newcomp) ) ] * 
                                                               rep(rep(sagw[2:8]/ agr, each=2 ),2) ) #CIN2 and 3 qalys women only
                                                      + sum(q[3,y] * x[c(fag(4, newcomp), fag(5, newcomp) )]* rep(rep(sagw[2:8]/ agr, each=2 ),2) ) #localised cancer
                                                      + sum(q[4,y] * x[c(fag(6, newcomp), fag(7, newcomp) )]* rep(rep(sagw[2:8]/ agr, each=2 ),2) ) #regional
                                                      + sum(q[5,y] * x[c(fag(8, newcomp), fag(9, newcomp) )]* rep(rep(sagw[2:8]/ agr, each=2 ),2) ) #distant
                                                      # - ( sum( (rme*Qme + rse*Qse)* (as.vector(pv[,1,vs])*g[length(agegroups)]* as.vector(N[,,length(agegroups) ])*yp ) ) ) #QALY of AE
                                                      ) ) 
    
    # cqoutcome[,8,i] <-  apply(cqoutcome[,1:5,i], 1, function(x) sum(x )  ) #with AE costs
    cqoutcome[,8,i] <-  apply(cqoutcome[,1:4,i], 1, function(x) sum(x )  ) #without AE costs
    
  }
  
  
  return(cqoutcome)
}

############### Multivariate sensitivity analysis on CE plane

#1. function that gives IC and IQ for each 4 strategies
# It includes the possibility to add 5 additionnal years with no vaccination health effect (but includes vaccination costs): SENS A.

CEA <- array(0, dim=c( length(vaccstrat),2, 2), 
             dimnames = list(vaccstrat , c("compared_baseline", "compared_no_vacc"), c("incremental_costs", "incremental_QALY") )  )

CEA_f <- function(dfinal2,nbsens, C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q){
  tempd <- fcostsvaccstrat_varyp2(dfinal2,nbsens,C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q)
  tempd <- tempd * popmod*sum(sagw[2:8])
  
  oCosts_b <- sum(tempd[startt:endt,8,1]* (exp(-disc*startt:endt) ) ) # for comparison with current vaccination strategy as baseline (50% V4v)
  
  oQALY_b <- sum(tempd[startt:endt,7,1] * (exp(-disc*startt:endt) ) )
  
  oCosts_nv <- sum( sum(tempd[1,1:3,1])* (exp(-disc*startt:endt) ) )  #for comparison to no vaccination (starting values)
  
  oQALY_nv <- sum ( (tempd[1,7,1] )* (exp(-disc*startt:endt) ) )  
  
  iCosts <- apply(tempd[,8,], 2, function(x) sum(x[startt:endt]*(exp(-disc*startt:endt) )  ) )
  iQALY <- apply(tempd[,7,], 2, function(x) sum(x[startt:endt]* (exp(-disc*startt:endt) )  ) )
  
  CEA[,1,1] <- (iCosts - oCosts_b)
  CEA[,1,2] <-  (iQALY - oQALY_b)
  
  CEA[,2,1] <- (iCosts - oCosts_nv)
  CEA[,2,2] <-  (iQALY - oQALY_nv)
  return(CEA) }

# New version (03.04.2019) with comparaison with same coverage rates
CEA <- array(NA, dim=c( length(vaccstrat),3, 2), 
             dimnames = list(vaccstrat , c("compared_56", "compared_80", "compared_no_vacc"), c("incremental_costs", "incremental_QALY") )  )

CEA_f <- function(dfinal2,nbsens, C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q){
  tempd <- fcostsvaccstrat_varyp2(dfinal2,nbsens,C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q)
  tempd <- tempd * popmod*sum(sagw[2:8])
  
  oCosts_56 <- sum(tempd[startt:endt,8,1]* (exp(-disc*startt:endt) ) ) # for comparison with current vaccination strategy as baseline (56% V4v)
  oCosts_80 <- sum(tempd[startt:endt,8,2]* (exp(-disc*startt:endt) ) ) # for comparison with current vaccination strategy as baseline (8% V4v)
  
  oQALY_56 <- sum(tempd[startt:endt,7,1] * (exp(-disc*startt:endt) ) )
  oQALY_80 <- sum(tempd[startt:endt,7,2] * (exp(-disc*startt:endt) ) )
  
  oCosts_nv <- sum( sum(tempd[1,1:3,1])* (exp(-disc*startt:endt) ) )  #for comparison to no vaccination (starting values), 
  oQALY_nv <- sum ( (tempd[1,7,1]+ tempd[1,6,1])* (exp(-disc*startt:endt) ) ) #total 
  
  iCosts <- apply(tempd[,8,], 2, function(x) sum(x[startt:endt]*(exp(-disc*startt:endt) )  ) )
  iQALY <- apply(tempd[,7,], 2, function(x) sum(x[startt:endt]* (exp(-disc*startt:endt) )  ) )
  
  CEA[3,1,1] <- (iCosts[3] - oCosts_56)
  CEA[3,1,2] <-  (iQALY[3] - oQALY_56)
  
  CEA[4,2,1] <- (iCosts[4] - oCosts_80)
  CEA[4,2,2] <-  (iQALY[4] - oQALY_80)
  
  CEA[,3,1] <- (iCosts - oCosts_nv)
  CEA[,3,2] <-  (iQALY - oQALY_nv)
  return(CEA) }


#3. take lsens number of different combinations of parameters.
# varied parameters are: "discounting", "treatement costs",  "vaccine costs", "utilities weights"
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return( c(alpha = alpha, beta = beta))
}

########### UNIVARIATE SENS ANALYSIS

ICER_fix <- array(0, dim=c( length(vaccstrat), 3), 
                  dimnames = list(vaccstrat ,  c("ICER_compared_56",  "ICER_compared_80","ICER_compared_to_novacc"  ) ))

icerf <- function(dfinal2,nbsens, C4v, C9v, disc,dis_c_sens, util_sens ){
  # disc <- discount_rate
  # dis_c_sens <- disease_cost
  # util_sens <- utilities
  tempd <- fcostsvaccstrat_varyp2(dfinal2,nbsens,C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q)
  
  tempd <- tempd * totpop*sum(sagw[2:8])
  
  oCosts_56 <- sum(tempd[startt:endt,8,1]* (exp(-disc*startt:endt) ) ) # for comparison with current vaccination strategy as baseline (50% V4v)
  oQALY_56 <- sum(tempd[startt:endt,7,1] * (exp(-disc*startt:endt) ) )
  
  oCosts_80 <- sum(tempd[startt:endt,8,2]* (exp(-disc*startt:endt) ) ) # for comparison with current vaccination strategy as baseline (50% V4v)
  oQALY_80 <- sum(tempd[startt:endt,7,2] * (exp(-disc*startt:endt) ) )
  
  oCosts_nv <- sum( sum(tempd[startt,1:3,1])* (exp(-disc*startt:endt) ) )  #for comparison to no vaccination (starting values)
  oQALY_nv <- sum ( (tempd[startt,7,1] )* (exp(-disc*startt:endt) ) ) 
  
  iCosts <- apply(tempd[,8,], 2, function(x) sum(x[startt:endt]*(exp(-disc*startt:endt) )  ) )
  iQALY <- apply(tempd[,7,], 2, function(x) sum(x[startt:endt]* (exp(-disc*startt:endt) )  ) )
  
  ICER_fix[,1] <- (iCosts - oCosts_56) / (iQALY - oQALY_56)
  ICER_fix[,2] <- (iCosts - oCosts_80) / (iQALY - oQALY_80)
  ICER_fix[,3] <- (iCosts - oCosts_nv) / (iQALY - oQALY_nv)
  
  return(ICER_fix)
}

