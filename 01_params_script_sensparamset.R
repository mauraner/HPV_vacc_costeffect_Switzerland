# CE parameter sets for natural history sensitivity analysis
# MR, 19.12.2018, updated 23.05.19
#directory, where to load parameter data from
dirl <- "tab_params/" 

# Make a defined number of parameter sets which have to be used for all 13 HPV types

#parameters that are varied are :
# HPV infection clearance rate in men and women (gamma), immunity decline rate (omega), progression from infection to cin2 (pri1) and to cin3 (pri2), 
# progression from cin2 to cin3 (pr23), progression of local cancer to regional (prlcrc) and from regional to distant (prrcdc), 
# probability of diagnostic (z), progrability of cure with treatement (pi) , mortality of undiagnosed (um) and diagnosed cancer (dm). 
# vaccine efficacy (ve), screening rates , screening sensitivity and screening specificity.
psn <- 2000

ptv <- c("alpha", "gamma_f", "gamma_m",  "omega", "tau", "pri1", "pri2", "pr23" , "pr2s", "pr3s", "prlcrc", "prrcdc", 
         "z1" ,"z2", "z3" ,"pi1", "pi2","pi3", "um1","um2","um3","dm1","dm2","dm3",
         "screening1", "screening2", "screening3", "screening4", "screening5","screening6", "screening7",
         "sensitivity", "specificity", "ve")

paramsets <- array(0, dim=c(psn, length(ptv), length(hpvtype)), 
                   dimnames = list( c(1:psn),ptv,  hpvtype)  ) 

## requires baseline parameters 
source("01_params_script_baselineparam.R")

####################################################### Retrieving parameter uncertainties

# for alpha, I take the CI given in the Durham paper
alpha_min <- c(0.48,0.27,0.46,0.24,0.55,0.33,0.23,0.18,0.086,0.36,0,0, mean(c(0.48,0.27,0.46,0.24,0.55,0.33,0.23,0.18,0.086,0.36)))
alpha_max <- c(0.65,0.57,0.67,0.59,0.97,0.63,0.56,0.95,0.26,0.77,0,0, mean(c(0.65,0.57,0.67,0.59,0.97,0.63,0.56,0.95,0.26,0.77)))


#for gamma in women, I take the IQR range, indicated in the paper (Jaisamaran 2013) (in months, to get rate per year: 1/ (gamma /12)
gamma_f_months <- c(17.11, 11.84, 13.8, 12.00, 11.77,11.77, 11.48, 11.77, 11.77 ,11.77, 8.26, 8.26, 11.77 )
gamma_f_min <- c(7.8,   6.2,   6.43, 6.20,  6.20,6.20,6.20, 6.20, 6.20 ,6.20, 5.97, 5.97, 6.20) 
gamma_f_max <- c(30.26, 23.11, 28.89, 21.90, 20.03,20.03, 23.31, 20.03, 20.03 ,20.03, 17.57, 17.57, 20.03 ) 

#for gamma men, I use IQR of exponential distribution based on the reported median estimates for min and max values of uniform distribution
gamma_m_min <- c()
gamma_m_max <- c()
for(i in 1:length(hpvtype)){
  t1 <- summary( (rexp(100000, (log(2)/ (gamma_a2[i]) ) ) ) )
  gamma_m_min[i]<-  t1[2]
  gamma_m_max[i]<- t1[5]
}

# for omega (waning immunity) I take the 95% PI from the paper (Johnson 2012)
omega_min <-c(0.011,0.010,0.010,0.010,0.010,0.011,0.010,0.010,0.013,0.010,0.10,0.10,0.10) #hpv6,11 and ohr unknown was set to min (mean 95%PI, posterior interval)
omega_max <-c(0.032,0.125,0.101,0.102,0.099,0.128,0.073,0.053,0.042,0.078,0.128,0.128,0.128)#hpv6,11 and ohr unknown was set to max

#relative risk of re-infection following clearance, from Castellsagué et al. 2014, for cervical cancer 
# I take the indicated CI from the paper for the indicated types (HPV 16 and 18), for the other, I take the min and max from these two combined
tau_c_min <- c(0.53, 0.75, rep(0.53,11)) 
tau_c_max <- c(0.78,1.19, rep(1.19,11) )

#for progression rate to cin2 and cin3, I take the CI of the multivariate analysis from the Jaisanarn 2013 paper (not indicated types are considered as other HR)
# proportion after 2y (turn into rate by -log(1-pri1)/2 )
pri1_min <- (42/4824) *  c(6.84, 2.4, 3.56, 6.34, rep(1.94,2), 2.08,rep(1.94,3),0,0,rep(1.94,1) ) 
pri1_max <-  (42/4824) *  c(12.51, 5.27  ,7.29, 13.18 , rep(3.57,2), 6.4, rep(3.57,3),0,0, rep(3.57,1) )

pri2_min <-  (7/4824) *  c(9.97, 1.78, 3.46, 9.45, rep(1.62,2), 1.17,rep(1.62,3),0,0,rep(1.62,1) ) 
pri2_max <- (7/4824) *  c(43.95, 12.58  ,17.63, 44.35 , rep(7.59,2), 16.97, rep(7.59,3),0,0, rep(7.59,1) ) 

#progression from CIN2 to CIN3: 95% CI of proportion of women showing progression after 3 years, 15% (9%-26%) tranform in rate by: -log(1-pr2c)/3
pr23_min <- 0.09
pr23_max <- 0.26

#regression from CIN2 to S and from CIN3 to S
# for pr2s, I take the 95%CI of the % cleared after 3 years, depending if it is HPV16/18 or other HPV (proportion after 3y, (turn into rate by -log(1-pr2s)/3 ) )
pr2s_min <- c( rep( 0.39  ,2), rep( 0.65 ,11) )
pr2s_max <- c( rep( 0.72  ,2), rep( 0.89 ,11) )
# for pr3s, I take the 95%CI of the percentage after 10y in the subset with persitent disease (proportion after 10y, (turn into rate by -log(1-pr3s)/10 ) )
pr3s_min <- (1-0.217 )
pr3s_max <- (1-0.434) 


#progression localised to regional and regional to distant, from Campos et al (2014), probabilities per months,
#to transform into rate: -log(1- p )/ (1/12), I vary by a factor two
prlcrc_min <- 0.020/2 #
prlcrc_max <- 0.020*2 #

prrcdc_min <- 0.025/2
prrcdc_max <- 0.025*2


### #stage specific rate of diagn., from Campos (2014) Web appendix, table 1, these are probabilities per months
#there is no CI around the parameters, which seemed to be based on a personal communication from another paper
# I vary by a factor two the  initial parameter c(0.0174,0.0735,0.1746), to transform into rate: -log(1- z )/ (1/12)
z_min <- c(0.0174,0.0735,0.1746)/2
z_max <- c(0.0174,0.0735,0.1746)*2


#stage specific probability of cure with treatement, from Durham who cites: Elbasha EH et al. 2007, values not found from the indicated reference, no CI
# I take +/- 10% (max of 1)
pi_min <- pi[,1]- 0.1*pi[,1]
pi_max <- pi[,1]+ 0.1*pi[,1]
pi_max[1] <- 1

#mortality of undiagnosed and diagnosed cancers (stage specific), from Campos et al. (2014). Like the rate of diagnosis, it 
# is given in probabilities per time (months). To transform into rate: -log(1- um )/ (1/12). I will vary the probabilities values by a factor 2

um_min <- c(0.0016,0.0095,0.0293)/2
um_max <- c(0.0016,0.0095,0.0293)*2

dm_min <- c(0.0009,0.0036,0.0076)/2
dm_max <- c(0.0009,0.0036,0.0076)*2

# For the screening, I take the 95% CI of binomial test based on the sample sizes per age groups
s1  <- read.csv( paste(dirl, "sceeningrates_2012SHS_v2.csv", sep="") )
s1 <- as.data.frame(s1[,2:13])
s2 <- array(0, dim=c( length(agegroups) ),  dimnames = list( agegroups) )
#
#rate is (-log(1-s2_min) )/1
s2_min <-  c( 0,0 , s1[1:5,11] )
s2_max<-  c( 0,0, s1[1:5,12] )


# For sensitivity and specificity of screening, I took the 95%CI from Bigras 2005
ssens_min  <- 0.486 #from Bigras 2005, CI48.6-68.2
ssens_max  <- 0.682 #from Bigras 2005, CI48.6-68.2
sspec_min <- 0.966 #from Bigras 2005, CI 96.6-97.2
sspec_max <- 0.972 #from Bigras 2005, CI 96.6-97.2


#Vaccine efficacy #hpv16:0.91-0.96, 18:0.91-0.97): from Lu et al 2011. 
#for the other hpv types: 0.967 CI 0.809-0.998) from Joura et al 2015
#ve <- c(0.94, 0.95, 0.967, 0.967, 0, 0,0.967,0,0.967, 0.967, 0.967, 0.967, 0)

ve_min <- c(0.91,0.91,rep(0.809,11))
ve_min[c(5,6,8,13)] <- 0 #non vaccine types
ve_max <- c(0.96,0.97,rep(0.998,11) )
ve_max[c(5,6,8,13)] <- 0

# --------------------------------------------------------------------------------------------------------------------------------
# create the parameter array, with sampling hypercube method, based on uniform distribution of max and min values,

#Based on script from Anthony Hauser, 11.12.2018 (sampling_hypercube_anthony.R)

require(lhs) #add the lhs library

# number of points (already defined above):
psn

#number of parameters to vary (partly defined above, but some vary by gender (gamma) or by age groups (screening)):
# nptvt <- length(ptv) + (length(sex)-1) + (length(agegroups)-1) + 4*(3-1) # z,pi, um and dm are for three cancer stages
ptv

set.seed(6242016)

# lhs<-randomLHS(psn,nptvt) #simulate
lhs<-maximinLHS(psn, length(ptv) ) #simulate

#first those who do not vary by HPV type 


paramsets[ , which(ptv=="pr23"),] <- -log(1-  (lhs[,which(ptv=="pr23")]*(pr23_max - pr23_min) + pr23_min) ) /3
paramsets[ , which(ptv=="prlcrc"),] <- -log(1- (lhs[,which(ptv=="prlcrc")]*(prlcrc_max - prlcrc_min) + prlcrc_min) ) / (1/12)
paramsets[ , which(ptv=="prrcdc"),] <- -log(1- (lhs[,which(ptv=="prrcdc")]*(prrcdc_max - prrcdc_min) + prrcdc_min) ) / (1/12)

for(i in 1:3){
  paramsets[ , which(ptv==paste("z",i, sep="")),] <- (-log(1- (lhs[,which(ptv==paste("z",i, sep=""))]*(z_max[i] - z_min[i]) + z_min[i]) ) ) / (1/12)
  paramsets[ , which(ptv==paste("pi",i, sep="")),] <- lhs[,which(ptv==paste("pi",i, sep=""))]*(pi_max[i] - pi_min[i]) + pi_min[i]
  paramsets[ , which(ptv==paste("um",i, sep="")),] <- (-log(1- (lhs[, which(ptv==paste("um",i, sep=""))]*(um_max[i] - um_min[i]) + um_min[i]) )) / (1/12)
  paramsets[ , which(ptv==paste("dm",i, sep="")),] <- (-log(1- (lhs[,which(ptv==paste("dm",i, sep=""))]*(dm_max[i] - dm_min[i]) + dm_min[i]) )) / (1/12)
}

for(j in 1:length(agegroups) ){
  paramsets[ , which(ptv==paste("screening",j,sep="")),] <-(-log(1- (lhs[,which(ptv==paste("screening",j,sep=""))]*(s2_max[j] - s2_min[j]) + s2_min[j]) )) /1
}

paramsets[ , which(ptv=="sensitivity"),] <- lhs[,which(ptv=="sensitivity")]*(ssens_max - ssens_min) + ssens_min #
paramsets[ , which(ptv=="specificity"),] <- lhs[,which(ptv=="specificity")]*(sspec_max - sspec_min) + sspec_min



for(i in 1:length(hpvtype)){
  paramsets[ , which(ptv=="alpha"),i]  <- lhs[,which(ptv=="alpha")]*(alpha_max[i]-alpha_min[i]) + alpha_min[i]
  paramsets[ , which(ptv=="gamma_f"),i] <-  1/( (lhs[,which(ptv=="gamma_f")]*(gamma_f_max[i] - gamma_f_min[i]) + gamma_f_min[i])/12)
  paramsets[ , which(ptv=="gamma_m"),i] <- 1/( (lhs[,which(ptv=="gamma_m")]*(gamma_m_max[i] - gamma_m_min[i]) + gamma_m_min[i])/12)
  paramsets[ , which(ptv=="omega"),i]  <- lhs[,which(ptv=="omega")]*(omega_max[i] - omega_min[i]) + omega_min[i]
  paramsets[ , which(ptv=="ve"),i] <-    lhs[,which(ptv=="ve")]*(ve_max[i] - ve_min[i]) + ve_min[i]
  paramsets[ , which(ptv=="tau"),i] <-    lhs[,which(ptv=="tau")]*(tau_c_max[i] - tau_c_min[i]) + tau_c_min[i]
  paramsets[ , which(ptv=="pri1"),i] <-  -log( 1- (lhs[,which(ptv=="pri1")]*(pri1_max[i] - pri1_min[i]) + pri1_min[i]) ) /2
  paramsets[ , which(ptv=="pri2"),i] <-  -log( 1- (lhs[,which(ptv=="pri2")]*(pri2_max[i] - pri2_min[i]) + pri2_min[i]) ) /2
  paramsets[ , which(ptv=="pr2s"),i] <-  -log( 1- (lhs[,which(ptv=="pr2s")]*(pr2s_max[i] - pr2s_min[i]) + pr2s_min[i]) ) /3
  paramsets[ , which(ptv=="pr3s"),i] <-  -log( 1- (lhs[,which(ptv=="pr3s")]*(pr3s_max - pr3s_min) + pr3s_min) ) /10
  
}  



hist(qbeta(lhs[,4],2,2) *(gamma_f_max[i] - gamma_f_min[i]) + gamma_f_min[i] )
# hist( qbeta(lhs[,3],2,2)*(alpha_max[i]-alpha_min[i]) + alpha_min[i])
summary( 1/ ( (lhs[,which(ptv=="gamma_f")]*(gamma_f_max[i] - gamma_f_min[i]) + gamma_f_min[i])/12 ) )
summary(  1/( (lhs[,which(ptv=="gamma_m")]*(gamma_m_max[i] - gamma_m_min[i]) + gamma_m_min[i])/12) )


dim(paramsets)

summary(paramsets[,which(ptv=="pr2s"),1])
summary(paramsets[,which(ptv=="pr3s"),1])

summary(paramsets[,which(ptv=="screening3"),1])

paramsets[1:10,which(ptv=="pr2s"),]

# summary(qbeta(lhs[,4],1,2.5) *(gamma_f_max[i] - gamma_f_min[i]) + gamma_f_min[i]) #normal distribution
# hist(qbeta(lhs[,4],1,2.5) *(gamma_f_max[i] - gamma_f_min[i]) + gamma_f_min[i]) #normal distribution

dirl <- "tab_params/" #where to load parameter data from

date_save <- gsub("-", "_", Sys.Date()) 


save(paramsets, file=paste(dirl,  "paramsets" , date_save, ".RData", sep="") )

dim(paramsets)
paramsets[1,2,1]


summary(paramsets[,1,1])


summary(paramsets[,2, ]  )

for(j in 1:length(hpvtype) ){
  for(i in 1:length(ptv)){
    print(ptv[i])
    print(summary(paramsets[,i,j]) )
    hist(paramsets[,i,j], main=paste(ptv[i], hpvtype[j], sep=": ") )
  }}

for(i in 1:length(agegroups)){
  print(agegroups[i])
  print(summary(paramsets[,5,1,i,]) )
  
}

summary(paramsets[,2,2,1,])

summary(paramsets[,4,1,1,])

######
# pre and posterior distribution of parameter values

for(i in 1:length(ptv)){
  print(ptv[i])
  print(round(summary(paramsets[,i,1])[c(4,1,6)],3))
  # print(summary(paramsets_ok[,i,1]))
}

dim(paramsets_ok)

for(i in 1:length(hpvtype)){
  print(hpvtype[i])
  print( round( summary(beta_ok[i,])[c(2,4,5)] ,2) )
  print(round(summary(pr3c_ok[i,])[c(2,4,5)] ,3) )
  print(round(summary(pr3c_ok[i,])[c(2,4,5)]*0.2 ,3) )
}

