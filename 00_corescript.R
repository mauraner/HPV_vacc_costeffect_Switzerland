#
#
# Core script of the HPV cost-effectiveness project in Switzerland


# Manuscript title: Impact and cost-effectiveness of nonavalent HPV vaccination in Switzerland: 
# insights from a dynamic transmission model
# Authors: Maurane Riesen, Johannes A. Bogaards, Nicola Low, Christian L. Althaus

# Maurane Riesen, final version May 2019
# ---------------------------------------------------------------------------------------------------

# The overall script is composed by a core script which is linked to 
# - R files to retrieve required parameters and a folder containing parameter tables (01_params_script_*, tab_params)
# - R files meant to be run on the high performance computing cluster UBELIX, used to do the 
#   calibration and sensitivity analysis, in two distinct steps (02_UBE_part*)
# - R files meant to be run on the high performance computing cluster UBELIX (02_UBE_part*)
# - R files containing the different functions (03_function_*)

# ------------------------------------------------------------------------------------------------------

########
# STEP 1: Prepare parameter sets
########
source("01_params_script_baselineparam.R")

# 01_params_* scripts and mainly 01_params_script_baselineparam.R and 01_params_script_sensparams.R
# are used to create 2000 different parameter sets by varying 34 parameters linked with the
# the natural history of cervical cancer, for all different types. The output parameterset is saved as:

# paramsets2019_01_26.RData

# Varied parameters:
# c("alpha", "gamma_f", "gamma_m",  "omega", "tau", "pri1", "pri2", "pr23" , "pr2s", "pr3s", "prlcrc", "prrcdc", 
#   "z1" ,"z2", "z3" ,"pi1", "pi2","pi3", "um1","um2","um3","dm1","dm2","dm3",
#   "screening1", "screening2", "screening3", "screening4", "screening5","screening6", "screening7",
#   "sensitivity", "specificity", "ve")

########
# STEP 2: Run the calibration, 1st step of the sensitivity analysis
########

#We calibrate the transmission probablity per partnership (beta) to our given HPV prevalence in 20-24 y.o women
# and the cancer progression rate from CIN3 to regional cervical cancer (pr3c), incl. pr2c which is assumed to be 0.2* pr3c

# This is done by first checking if the parameter set is able to reach the targetted HPV prevalence 
# (using beta = 1 and assuring that the prevalence is higher than the targetted prevalence)
# retrieving the correct beta for each hpv type by running the simulation on different beta values
# and using a linear regression between the two beta estimates achieving prevalences closest to the targetted prevalence

# The calibrated beta is then used to retrieve the pr3c using the same method, but by fitting the crude Swiss 
# cervical cancer incidence on the modelled age groups

# This process is coded in 02_UBE_sens1.R which is run on the UBELIX servers, using the parameter set paramsets2019_01_26.RData
# The output files of this first step of the sensitivity analysis are copied in the folder sens_step1


########
# STEP 3: Retrieve working parameter sets and run the second step of the sensitivity analysis
########  

# The simulation for each parameter sets (incl all HPV types) are stored in the folder sens_step1,
# now we retrieve the parameter sets which achieved calibration and use those to continue the analysis.
# This is done using script: 01_params_script_workingparamsets.R which checkt the output and stores the 
# calibrated values for beta and pr3c, and also makes the parameterset tables before and after calibration

# It creates a vector of calibrated beta and pr3c parameters and a new parameterset with only the working sets
# beta_ok<date>.RData, pr3c_ok<date>.RData and paramsets_ok<date>.RData which are used to run the script:
# 02_UBE_sens2.R on Ubelix (step2_sens.job). 

########
# STEP 4: Load the final simulations and get the results for health effects
########  
dirl <- "tab_params/"

load(paste(dirl, "paramsets_ok2019_05_28.RData", sep="") )
load(paste(dirl, "beta_ok2019_05_28.RData", sep="") )
load(paste(dirl, "pr3c_ok2019_05_28.RData", sep="") ) 

load(paste(dirl, "paramsets_ok2019_06_04.RData", sep="") )
load(paste(dirl, "beta_ok2019_06_04.RData", sep="") )
load(paste(dirl, "pr3c_ok2019_06_04.RData", sep="") ) 

#directories
dirl <- "tab_params/" #where to load parameter data from
plotdir <- "plots/" #where to save plots
loaddir <- "sens_step2/190603/p" #where the UBE data was stored 
# dirlts <- "O:/Cost_effectiveness project/Preparation_drafts/tab/" #where to save the latex tables to
dirlts <- "C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/" #where to save the latex tables to

source("01_params_script_baselineparam.R")
source("03_functions_costs.R")

ptv <- c("alpha", "gamma_f", "gamma_m",  "omega", "tau", "pri1", "pri2", "pr23" , "pr2s", "pr3s", "prlcrc", "prrcdc", 
         "z1" ,"z2", "z3" ,"pi1", "pi2","pi3", "um1","um2","um3","dm1","dm2","dm3",
         "screening1", "screening2", "screening3", "screening4", "screening5","screening6", "screening7",
         "sensitivity", "specificity", "ve")


# recall parameters used for the simulations
time <- seq(0, 200, 1)
diffv <- c(0.56, 0.8)
vaccstrat <-  c( paste( rep("V4v", each =length(diffv) ), diffv*100, sep="_"), paste( rep("V9v", each =length(diffv) ), diffv*100, sep="_") )
y <-  which(dis_outcome=="cervix")

# Load simulations

load( paste(loaddir,6,".RData", sep="")) # test

#indicate length of working parameter sets
nsens_ok <- 1116



dfinal3 <- array(NA, dim=c(nsens_ok, length(time),length(col_nam2),length(vaccstrat) ), dimnames = list(1:nsens_ok,time,col_nam2, vaccstrat) )  

for(i in 1:nsens_ok){
  if( file.exists(paste(loaddir, i,".RData", sep="") ) ) {
    load( paste(loaddir, i,".RData", sep="") )
    dfinal3[i,,,] <- tosave[[1]]
    next
  }}  

# ************ plot cancer, CIN and mortality over time *****************************
dim(dfinal3)

dfinal3[6,3,,1]

table(is.na(dfinal3[,2,2,1]))

which(is.na(dfinal3[,2,2,1]))

totmodinc <- array(0, dim=c(nsens_ok, length(time), length(vaccstrat), 3),  dimnames = list(1:nsens_ok, time, vaccstrat , c("incidence", "50%", "95%")  ) )
tempin <- c()
for(kn in 1:nsens_ok){
  for(j in 1:length(time) ){
    for( i in 1:length(vaccstrat)){
      for(k in 0:(length(agegroups)-1) ){
        tempin[k+1] <- sum(dfinal3[kn,j, 1+c(k*(10*4)+4,k*(10*4)+6 ,k*(10*4)+8, k*(10*4)+(10*2)+4,k*(10*4)+(10*2)+6 ,k*(10*4)+(10*2)+8) ,i]   
                           * rep( c(paramsets_ok[kn,which(ptv=="z1"),1],
                                    paramsets_ok[kn,which(ptv=="z2"),1],
                                    paramsets_ok[kn,which(ptv=="z3"),1]), 2) / (sr*agr[1+k])  * sagw[2+k])
      }
      totmodinc[kn,j,i,1] <- sum(tempin  )
    }} 
}

popmod <- totpop
totmodinc[1:20,1,1,1] * popmod*sr

denom <-   popmod*sr*sum(sagw[2:8])
totmodinc[100:120,1,1,1] * denom

summary(totmodinc[,1,1,1]*1e5)

which(totmodinc[,1,1,1]>7.200e-05)

totmodinc[1108,1,1,1]

# totmodinc[c(6,74,127,177,431,497,520,647,650,657,701,920,1068),,4,1] <- NA

# I take the scenario within 50% and 95% quantiles after 50y, excluding 0
crosst <- 50

for(i in 1:length(vaccstrat)){
  totmodinc[,,i,2] <- cut(totmodinc[,crosst,i,1], quantile(totmodinc[,crosst,i,1] ,probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE)   )
  totmodinc[,,i,3] <- cut(totmodinc[,crosst,i,1], quantile(totmodinc[,crosst,i,1] ,probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE)   )
}

table( totmodinc[,1,1,1]==0 )
table( totmodinc[,1,1,3])

table(totmodinc[,1,1,2])
min(totmodinc[,1,4,1] * denom, na.rm=TRUE)

new_totmodinc <- array(0, dim=c(5, length(time), length(vaccstrat)),  dimnames = list( c("min50", "max50" ,"min95", "max95","med"), time, vaccstrat   ) )

for(j in 1:length(time)){
  for(i in 1:length(vaccstrat)){
    
    new_totmodinc[1,j,i] <-  min(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,2]) ] )
    new_totmodinc[2,j,i] <-  max(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,2]) ] ) 
    
    new_totmodinc[3,j,i] <-  min(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,3])  ] )
    new_totmodinc[4,j,i] <-  max(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,3])   ]) 
    
    new_totmodinc[5,j,i] <-  median(totmodinc[,j,i,1], na.rm=TRUE) 
    
  }}


new_totmodinc[,80,]*denom
new_totmodinc[,1,]*denom

######### plotting:
# col3 <- colorRampPalette(c("darkblue", "lightblue"))(length(diffv))
# col4 <- colorRampPalette(c("darkgreen", "lightgreen"))(length(diffv))
# col1 <- c( col3, col4)
al1 <- 0.6
al2 <- 1  
col1 <- c(rgb(0.44,0.55,0.86, alpha=al1, maxColorValue = 1),  rgb(0.27,0.31,0.43, alpha=al1, maxColorValue = 1) ,  
          rgb(0.4,0.74,0.36, alpha=al1, maxColorValue = 1), rgb(0.27,0.42,0.25, alpha=al1, maxColorValue = 1)   )

col2 <- c(rgb(0.44,0.55,0.86, alpha=al2, maxColorValue = 1),   rgb(0.27,0.31,0.43, alpha=al2, maxColorValue = 1)  ,  
          rgb(0.4,0.74,0.36, alpha=al2, maxColorValue = 1), rgb(0.27,0.42,0.25, alpha=al2, maxColorValue = 1)   )

plott <- 100
stt <- 0
denom <- 1e5
denom <- popmod*sr*sum(sagw[2:8])

pdf(paste(plotdir, "CCincidence_sens", denom,"_", Sys.Date(), ".pdf", sep=""), width=8, height=6)
par(mar=c(5,7,7,2) ) 
rf <- 50

plot(NA, xlim=c(stt,plott ), ylim=c(0,round(max(totmodinc[,,,1]*denom, na.rm=TRUE)/rf)*rf ),
     xlab= "years",
     ylab="incident \ncervical cancer cases",
     cex.lab=2, frame.plot=FALSE)

for(i in c(1,2,3,4)){
  polygon( c(stt:plott, rev(stt:plott) ), c(new_totmodinc[2,(stt+1):(plott+1),i]*denom, rev(new_totmodinc[1,(stt+1):(plott+1),i]*denom) )
           ,col = col1[i], border = NA)
  
  polygon( c(stt:plott, rev(stt:plott) ), c(new_totmodinc[4,(stt+1):(plott+1),i]*denom, rev(new_totmodinc[3,(stt+1):(plott+1),i]*denom) )
           ,col = col1[i], border = NA)
  lines(stt:plott,new_totmodinc[5,(stt+1):(plott+1),i]*denom, col= col2[i], lwd=2 )
}

legend(70,round(max(totmodinc[,,,1]*denom, na.rm=TRUE)/rf)*rf, legend=c( paste("V4v:", (diffv*100), "%", sep="") , paste("V9v:", (diffv*100), "%", sep="")), 
       lty=1 , col=col1, lwd=rep(c(15,7,1), each=4), bty="n", title="Vaccine type\nand coverage" )

legend(70,round(max(totmodinc[,,,1]*denom, na.rm=TRUE)/rf)*rf, legend=c( "","","",""), 
       lty=1 , col=col1, lwd=10, bty="n", title="\n" )
legend(70,round(max(totmodinc[,,,1]*denom,na.rm=TRUE)/rf)*rf, legend=c( "","","",""), 
       lty=1 , col=col2, lwd=1, bty="n", title="\n" )

dev.off()


#decrease after 100 years *************************************
timt <- 101
new_totmodinc[5,timt,]

round((new_totmodinc[5,1,] * sum(sagw[2:8]))*1e5,1)
round((new_totmodinc[c(1:2,5),timt,] * sum(sagw[2:8]))*1e5,1)

round(1- (new_totmodinc[5,timt,]/new_totmodinc[5,1,]),2)
round((1- (new_totmodinc[3:5,timt,]/new_totmodinc[3:5,1,]) )*100,1)
round((1- (new_totmodinc[c(1:2,5),timt,]/new_totmodinc[c(1:2,5),1,]) )*100,1)

# Number of averted cases over 50y *************************************
dim(dfinal2)

outcomoi <- c("CIN diagnosed cases", "cervical cancer cases",  "cervical cancer deaths")
forplot <- c("min50", "max50" ,"min95", "max95","med")


for(i in 1:length(vaccstrat)){
  totmodinc[,,i,2] <- cut(totmodinc[,crosst,i,1], quantile(totmodinc[,crosst,i,1] ,probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE)   )
  totmodinc[,,i,3] <- cut(totmodinc[,crosst,i,1], quantile(totmodinc[,crosst,i,1] ,probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE)   )
}

table( totmodinc[,1,1,1]==0 )
table( totmodinc[,1,1,2])

summary(totmodinc[,1,1,2])
min(totmodinc[,1,4,1] * denom)

new_totmodinc <- array(0, dim=c(5, length(time), length(vaccstrat)),  dimnames = list( c("min50", "max50" ,"min95", "max95","med"), time, vaccstrat   ) )

for(j in 1:length(time)){
  for(i in 1:length(vaccstrat)){
    
    new_totmodinc[1,j,i] <-  min(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,2]) ] )
    new_totmodinc[2,j,i] <-  max(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,2]) ] ) 
    
    new_totmodinc[3,j,i] <-  min(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,3])  ] )
    new_totmodinc[4,j,i] <-  max(totmodinc[,j,i,1][!is.na(totmodinc[,crosst,i,3])   ]) 
    
    new_totmodinc[5,j,i] <-  median(totmodinc[,j,i,1]) 
    
  }}
nsim
simn <- nsens_ok

crosst <- 50

totcasesinSwitz <- array(0, dim=c(length(time), length(vaccstrat), length(outcomoi), simn),  
                         dimnames = list(time, vaccstrat ,outcomoi, 1:simn)  )
totcas_plot <- array(0, dim=c(length(time), length(vaccstrat), length(outcomoi), length(forplot)),  
                     dimnames = list(time, vaccstrat , outcomoi, forplot)  )


#if 10-14 agegroup included:

for(j in 1:length(vaccstrat)){
  for(nsim in 1:simn ){
    for(tt in 1:length(time)){
      cintreatpt <- rep(0,length(hpvtype))
      cincpt <- rep(0,length(hpvtype))
      dccincpt <- rep(0,length(hpvtype))
      totcasesinSwitz[tt,j,1,nsim] <- sum( rep(rep(c(paramsets_ok[nsim,which(ptv=="screening1"),1],
                                                     paramsets_ok[nsim,which(ptv=="screening2"),1],
                                                     paramsets_ok[nsim,which(ptv=="screening3"),1],
                                                     paramsets_ok[nsim,which(ptv=="screening4"),1],
                                                     paramsets_ok[nsim,which(ptv=="screening5"),1],
                                                     paramsets_ok[nsim,which(ptv=="screening6"),1],
                                                     paramsets_ok[nsim,which(ptv=="screening7"),1] #only if 10-14 are included
      ), each=length(riskgroups)),2)
      *paramsets_ok[nsim,which(ptv=="sensitivity"),1]* (dfinal3[nsim,tt, c(fag(2,newcomp),fag(3,newcomp)) ,j ] / 
                                                          (sr* rep( rep(agr, each=length(riskgroups) ),2 ) ) )
      * rep( rep(sagw[2:8], each=length(riskgroups)  ),2)  ) 
      
      totcasesinSwitz[tt,j,2,nsim] <-  sum( rep( c(paramsets_ok[nsim,which(ptv=="z1"),1],
                                                   paramsets_ok[nsim,which(ptv=="z2"),1],
                                                   paramsets_ok[nsim,which(ptv=="z3"),1]), each=length(agegroups)*length(riskgroups))   
                                            * (dfinal3[nsim,tt, c(fag(4,newcomp),fag(6,newcomp),fag(8,newcomp)) , j]/
                                                 (sr*rep(rep(agr,each=length(riskgroups) ),3) ) ) * rep( rep( sagw[2:8], each=length(riskgroups) ),3) )
      
      
      totcasesinSwitz[tt,j,3,nsim] <-  sum( c( rep( c(paramsets_ok[nsim,which(ptv=="um1"),1],
                                                      paramsets_ok[nsim,which(ptv=="um2"),1],
                                                      paramsets_ok[nsim,which(ptv=="um3"),1]), each=length(agegroups)*length(riskgroups) ), 
                                               rep( c(paramsets_ok[nsim,which(ptv=="dm1"),1],
                                                      paramsets_ok[nsim,which(ptv=="dm2"),1],
                                                      paramsets_ok[nsim,which(ptv=="dm3"),1]), each=length(agegroups)*length(riskgroups)  ) )
                                            * (dfinal3[nsim,tt, c(fag(4,newcomp),fag(6,newcomp),fag(8,newcomp), fag(5,newcomp),fag(7,newcomp),fag(9,newcomp)),j]/
                                                 (sr*rep(rep(agr,each=length(riskgroups) ),6 ) ) * rep( rep( sagw[2:8], each=length(riskgroups)),6 ) ) ) 
      
    }}
  
  for(k in 1:length(outcomoi)){
    totcas_plot[,j,k,1] <- apply(totcasesinSwitz[,j,k,],1, function(x) min(x[!is.na(cut(totcasesinSwitz[crosst,j,k,], quantile(totcasesinSwitz[crosst,j,k,] ,probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE) )  )] ))
    totcas_plot[,j,k,2] <- apply(totcasesinSwitz[,j,k,],1, function(x) max(x[!is.na(cut(totcasesinSwitz[crosst,j,k,], quantile(totcasesinSwitz[crosst,j,k,] ,probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE) )  )] ))
    
    totcas_plot[,j,k,3] <- apply(totcasesinSwitz[,j,k,],1, function(x) min(x[!is.na(cut(totcasesinSwitz[crosst,j,k,], quantile(totcasesinSwitz[crosst,j,k,] ,probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE) )  )] ))
    totcas_plot[,j,k,4] <- apply(totcasesinSwitz[,j,k,],1, function(x) max(x[!is.na(cut(totcasesinSwitz[crosst,j,k,], quantile(totcasesinSwitz[crosst,j,k,] ,probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE) )  )] ))
    
    totcas_plot[,j,k,5] <- apply(totcasesinSwitz[,j,k,],1,function(x)    median(x, na.rm=TRUE) )
  }
}

summary(totcasesinSwitz[1,1,1,])
summary(totcasesinSwitz[1,1,2,])
summary(totcasesinSwitz[1,1,3,])

summary(totcasesinSwitz[1,4,3,]*denom)

summary(totcasesinSwitz[1,4,2,]*denom)

summary(totcasesinSwitz[1,4,1,]*denom)

summary(totcasesinSwitz[1,4,3,]*denom)

# denom1 <- popmod * sr * 1e6

totcas_plot_2 <- totcas_plot
totcas_plot_2 <- totcas_plot_2 *denom

totcasesinSwitz_2  <-  totcasesinSwitz
totcasesinSwitz_2 <- totcasesinSwitz_2 *denom
summary(totcasesinSwitz_2[1,1,1,])

###
al1 <- 0.4
al2 <- 0.8

col1 <- c(rgb(0.44,0.55,0.86, alpha=al1, maxColorValue = 1),  rgb(0.27,0.31,0.43, alpha=al1, maxColorValue = 1) ,  
          rgb(0.4,0.74,0.36, alpha=al1, maxColorValue = 1), rgb(0.27,0.42,0.25, alpha=al1, maxColorValue = 1)   )

col2 <- c(rgb(0.44,0.55,0.86, alpha=al2, maxColorValue = 1),   rgb(0.27,0.31,0.43, alpha=al2, maxColorValue = 1)  ,  
          rgb(0.4,0.74,0.36, alpha=al2, maxColorValue = 1), rgb(0.27,0.42,0.25, alpha=al2, maxColorValue = 1)   )

col3 <- c(rgb(0.65,0.65,0.68, alpha=al1, maxColorValue = 1), rgb(0.44,0.44,0.46, alpha=al2, maxColorValue = 1) )

# pdf(paste(plotdir, "CCdeath_inc", Sys.Date(), ".pdf", sep=""), width=9, height=6)
# pdf(paste(plotdir, "CCcin_inc", Sys.Date(), ".pdf", sep=""), width=9, height=6)
# pdf(paste(plotdir, "CC_inc", Sys.Date(), ".pdf", sep=""), width=9, height=6)

pdf(paste(plotdir, "CCcancerMortalityandCINtreatements", Sys.Date(), ".pdf", sep=""), width=6, height=6)
# par(mfrow=c(3,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(oma = c(3, 3, 2, 1)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(2, 2, 1, 1)) # make the plots be closer together
rf <- c(1000,50, 50)

for(j in 1:3 ){
  j1 <- c(2,1,3)[j]
  plot(NA, xlim=c(0,length(1:plott) ), ylim=c(0,round(max(totcas_plot_2[,,j1,])/rf[j1] )*rf[j1] ) , 
       ylab= 'years after vaccination onset',
       xlab=  'number of incident cases per year',
       frame.plot=FALSE,  main=outcomoi[j1])
  
  for(i in 1:length(vaccstrat)){
    # polygon( c(stt:plott, rev(stt:plott) ), c(totcas_plot_2[(stt+1):(plott+1),i,j1,2], rev(totcas_plot_2[(stt+1):(plott+1),i,j1,1]) )
    # ,col = col1[i], border = NA)
    
    polygon( c(stt:plott, rev(stt:plott) ), c(totcas_plot_2[(stt+1):(plott+1),i,j1,4], rev(totcas_plot_2[(stt+1):(plott+1),i,j1,3]) )
             ,col = col1[i], border = NA)
  }
  for(i in 1:length(vaccstrat)){
    lines(stt:plott,totcas_plot_2[(stt+1):(plott+1),i,j1,5], col= col2[i], lwd=3 )
  }
  
  if(j==1){
    legend(75,round(max(totcas_plot_2[,,j1,])/rf[j1] )*rf[j1], legend=c( paste("V4v:", (diffv*100), "%", sep="") , paste("V9v:", (diffv*100), "%", sep="")), 
           pch=15 , col=col1, pt.cex=2 , bty="n", title="Vaccine type\nand coverage" )
    
    legend(0,30, legend=c( "" ), 
           lwd=30 , col=col3[1], seg.len=1.5 , bty="n" )
    # legend(-1.5,30, legend=c( "" ), 
    # lwd=15 , col=col3[1], seg.len=2.5 , bty="n" )
    legend(-2,30, legend=c( "" ), 
           lwd=2 , col=col3[2], seg.len=3.4, bty="n" )
    
    legend(8,32, legend=c("95% of the sim.", "median",
                          # "50% of the sim."
                          "" ), bty="n", cex=0.5 )
  }
}
mtext('years after vaccination onset', side = 1, outer = TRUE, line = 1, cex=1.5)
mtext('number of incident cases per year', side = 2, outer = TRUE, line = 1, cex=1.5)
dev.off()

#differences in number of cin treatements and deaths between strategies and baseline
# vaccstrat2 <- c(  "V4v: 0.5" ,"V4v: 0.6", "V4v: 0.7", "V9v: 0.5" ,"V9v: 0.6", "V9v: 0.7")
vaccstrat_lat <- c("$V4v_{56}$", "$V4v_{80}$" ,"$V9v_{56}$" ,"$V9v_{80}$")

casesaverted <- array(0, dim=c(length(vaccstrat), length(outcomoi)),  dimnames = list( vaccstrat_lat ,outcomoi)  )
casesaverted2 <- casesaverted

for( i in 2:length(vaccstrat)){
  temp1 <- summary(colSums(totcasesinSwitz_2[2:102,1,1,]) - colSums(totcasesinSwitz_2[2:102,i,1,]))
  temp2 <- summary(colSums(totcasesinSwitz_2[2:102,1,2,]) - colSums(totcasesinSwitz_2[2:102,i,2,]))
  temp3 <- summary(colSums(totcasesinSwitz_2[2:102,1,3,]) - colSums(totcasesinSwitz_2[2:102,i,3,]))
  casesaverted[i,1] <- paste( round( temp1[3],0 )," (IQR:", round( temp1[2],0 ), "-", round( temp1[5],0 ),")", sep="")
  casesaverted[i,2] <- paste( round( temp2[3],0 )," (IQR:", round( temp2[2],0 ), "-", round( temp2[5],0 ),")", sep="") 
  casesaverted[i,3] <- paste( round( temp3[3],0 )," (IQR:", round( temp3[2],0 ), "-", round( temp3[5],0 ),")", sep="") 
  
}
totcasesinSwitz_3 <- totcasesinSwitz_2

# totcasesinSwitz_3 <- apply(totcasesinSwitz_2[,1,1,],2, function(x) rep(x[1], 201) )

#########Re-do averted cases based on comparison with similar vaccination uptakes

casesaverted_samecov <- array(0, dim=c(2, length(outcomoi)),  
                              dimnames = list( vaccstrat_lat[3:4]  ,outcomoi)  )

for( i in 1:2){
  temp1 <- summary(colSums(totcasesinSwitz_2[2:102,i,1,]) - colSums(totcasesinSwitz_2[2:102,i+2,1,]))
  temp2 <- summary(colSums(totcasesinSwitz_2[2:102,i,2,]) - colSums(totcasesinSwitz_2[2:102,i+2,2,]))
  temp3 <- summary(colSums(totcasesinSwitz_2[2:102,i,3,]) - colSums(totcasesinSwitz_2[2:102,i+2,3,]))
  casesaverted_samecov[i,1] <- paste( round( temp1[3],0 )," (IQR:", round( temp1[2],0 ), "-", round( temp1[5],0 ),")", sep="")
  casesaverted_samecov[i,2] <- paste( round( temp2[3],0 )," (IQR:", round( temp2[2],0 ), "-", round( temp2[5],0 ),")", sep="") 
  casesaverted_samecov[i,3] <- paste( round( temp3[3],0 )," (IQR:", round( temp3[2],0 ), "-", round( temp3[5],0 ),")", sep="") 
}


for( i in 1:length(vaccstrat)){
  temp1 <- summary( colSums( apply(totcasesinSwitz_2[,1,1,],2, function(x) rep(x[1], 101) ) )  - colSums(totcasesinSwitz_2[2:102,i,1,]))
  temp2 <- summary(colSums( apply(totcasesinSwitz_2[,1,2,],2, function(x) rep(x[1], 101) ) )  - colSums(totcasesinSwitz_2[2:102,i,2,]))
  temp3 <- summary(colSums( apply(totcasesinSwitz_2[,1,3,],2, function(x) rep(x[1], 101) ) )  - colSums(totcasesinSwitz_2[2:102,i,3,]))
  
  casesaverted2[i,1] <- paste( round( temp1[3],0 )," (IQR:", round( temp1[2],0 ), "-", round( temp1[5],0 ),")", sep="")
  casesaverted2[i,2] <- paste( round( temp2[3],0 )," (IQR:", round( temp2[2],0 ), "-", round( temp2[5],0 ),")", sep="") 
  casesaverted2[i,3] <- paste( round( temp3[3],0 )," (IQR:", round( temp3[2],0 ), "-", round( temp3[5],0 ),")", sep="") 
}

casesaverted2

sum(rep(totcasesinSwitz[1,1,1], length(2:102)))

# write.table( casesaverted, paste("C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/cases_avertedall_1_", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(  casesaverted2, paste("C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/cases_avertedall_2_", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
write.table(  casesaverted_samecov, paste("C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/cases_averted_compcov", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)


#################################################### Cost/QALY SENS

date_save <- gsub("-", "_", Sys.Date()) 

save(dfinal3, file=paste( "dfinal3_", date_save, ".RData",sep="") )


########
# STEP 5: Get results for the cost-effectiveness analysis
########  

source("01_params_script_baselineparam.R")
source("03_functions_costs.R")
source("01_params_script_costsqaly_baselineparam.R")

#Load cost/qaly parameter set
load("tab_params/paramsets_CQok2019_06_04.RData")

# Load final simulation from sensitivity analysis and add Cost and Qaly sens analysis
load( paste("dfinal3_2019_06_11",".RData", sep="") )
dim(dfinal3)

ptv_cq <- c("disc", "Ct_1", "Ct_2",  "Ct_3", "Ct_4", "Ct_5", "Cpt", "Cfp" , "Cme", "Cse", "C4v", "C9v", 
            "q_1" ,"q_2", "q_3" ,"q_4", "q_5")

dirl <- "tab_params/"
plotdir <- "plots/" #where to save plots

# add costs aspect to model calculated without the costs
diffv <- c(0.56, 0.8)
vaccstrat <-  c( paste( rep("V4v", each =length(diffv) ), diffv*100, sep="_"), paste( rep("V9v", each =length(diffv) ), diffv*100, sep="_") )
time <- seq(0, 200, 1)

#give yp which is the weight to be applied on the incoming population in order to scale it to the Swiss population structure
popmod <- totpop # population size end of 2016

# in 2016, from the total population residing in Switzerland (all nationality, both sexes, permanently residing population)
# in the age group 11-14 years (which is the vaccination target age) the population sizes are 11y:81149, 12y:81083, 13y: 79781, 14y:81032, mean: 80761.25
yp <- 80761.25 / (g[length(agegroups)]* sum( as.vector(N[,,length(agegroups) ]) )  * popmod *sum(sagw[2:8]) )
g[length(agegroups)]* sum( as.vector(N[,,length(agegroups) ]) ) * yp  * popmod *sum(sagw[2:8])
#check with Swiss population age structure
(sagw * popmod)[1] / 10 

#calculates also costs and qaly, for overall population (weighted according to Swiss pop structure)
diffcosts <- c("pap_costs", "CIN2CIN3_costs", "cancer_costs", "vacc_costs", "AE_costs", "AE_QALY", "TOT_QALY", "TOT_costs")


######################################### ICER

############### Multivariate sensitivity analysis on CE plane,
# It randomly attributes a cost/qaly parameter to each scenario from the sensitivity analysis

#1. function that gives IC and IQ for each 4 strategies: CEA_f

#look, after 100y vaccination start, the ICER of different vaccination strategies (compared to baseline/current  situation)
# integrate over defined time period (area under curve)
#for the vaccination cohort life time (15-84y.o) i.e 69 years
endt <- 102
startt <-2

#discouting values
disc <- 0.03

#2. empty array to store results
#paramset for Cost and qaly
dim(paramsets_ok)
lsens <- nsens_ok

CEA_sens <- array(0, dim=c(  length(vaccstrat), 2, 2, lsens ), 
                  dimnames = list(vaccstrat , c("compared_baseline", "compared_no_vacc"), 
                                  c("incremental_costs", "incremental_QALY"),1:lsens))

CEA_sens <- array(NA, dim=c(  length(vaccstrat), 3, 2, lsens ), 
                  dimnames = list(vaccstrat , c("compared_56", "compared_80", "compared_no_vacc"),
                                  c("incremental_costs", "incremental_QALY"),1:lsens))

dim(paramsets_CQok)

pb <- txtProgressBar(min = 0, max = lsens, style = 3)
for(nbsens in 1:lsens){
  #discounting rate
  disc <-  paramsets_CQok[nbsens,which(ptv_cq=="disc")]
  #disease treatement costs
  Ct[1,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="Ct_1")]
  Ct[2,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="Ct_2")]
  Ct[3,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="Ct_3")]
  Ct[4,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="Ct_4")]
  Ct[5,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="Ct_5")]
  Cpt    <-  paramsets_CQok[nbsens,which(ptv_cq=="Cpt")]  #costs per negative cercival screen
  Cfp    <-  paramsets_CQok[nbsens,which(ptv_cq=="Cfp")] #costs of follow-up of false positive
  Cme    <-  paramsets_CQok[nbsens,which(ptv_cq=="Cme")] #cost per mild adverse effect (ae)
  Cse    <-  paramsets_CQok[nbsens,which(ptv_cq=="Cse")] #cost per severe adverse effect (ae)
  #vaccination costs and price diff 
  C4v <-   paramsets_CQok[nbsens,which(ptv_cq=="C4v")]
  C9v <- paramsets_CQok[nbsens,which(ptv_cq=="C9v")]
  # utilities weights
  q[1,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="q_1")]
  q[2,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="q_2")]
  q[3,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="q_3")]
  q[4,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="q_4")]
  q[5,1]  <-  paramsets_CQok[nbsens,which(ptv_cq=="q_5")]
  ##run function and save in array
  CEA_sens[,,,nbsens] <-  CEA_f(dfinal3[nbsens,,,],nbsens,C4v, C9v, Ct, Cpt, Cfp, Cme,Cse, q)
  setTxtProgressBar(pb, nbsens)
}


# save( CEA_sens , file="CEA_sensv2.RData")
# load("CEA_sensv2.RData")

################################## Plot
require(reshape2)
require(ggplot2)

#new version: compare with similar vaccination coverages (03.04.2019)
mCEA_sens56 <- melt(CEA_sens[,1,,])
mCEA_sens56 <- mCEA_sens56[-which(is.na(mCEA_sens56[,4])),]

mCEA_sens80 <- melt(CEA_sens[,2,,])
mCEA_sens80 <- mCEA_sens80[-which(is.na(mCEA_sens80[,4])),]

mCEA_sens1 <- rbind(mCEA_sens56, mCEA_sens80)

head(mCEA_sens1[600:1000,])
head(mCEA_sens80)

#1st for comparison with current situation

# mCEA_sens1 <- melt(CEA_sens[,1,,])
mCEA_sens1b <- dcast(mCEA_sens1, formula = Var3+Var1~ Var2 )
mCEA_sens1b$incremental_costs <- mCEA_sens1b$incremental_costs / 1e6 #costs in mio and qaly in tousands

#take 50% and 95% of the values for costs or qalys
for(i in 2:length(vaccstrat)){
  mCEA_sens1b$cut_cost_50[ mCEA_sens1b$Var1==vaccstrat[i] ] <- cut( mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],3 ],  quantile(mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],3 ], probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE) )
  mCEA_sens1b$cut_cost_95[ mCEA_sens1b$Var1==vaccstrat[i] ] <- cut( mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],3 ],  quantile(mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],3 ],  probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE) )
  mCEA_sens1b$cut_qaly_50[ mCEA_sens1b$Var1==vaccstrat[i] ] <- cut( mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],4 ],quantile(mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],4 ],probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE) ) 
  mCEA_sens1b$cut_qaly_95[ mCEA_sens1b$Var1==vaccstrat[i] ] <- cut( mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],4 ],quantile(mCEA_sens1b[ mCEA_sens1b$Var1==vaccstrat[i],4 ],probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE) ) 
}

mCEA_sens1b$cut_costandqaly_IQR <- ifelse(!is.na(mCEA_sens1b$cut_cost_50) & !is.na(mCEA_sens1b$cut_qaly_50), "50%", NA) 

#2cd for comparison with no vacc
# mCEA_sens2 <- melt(CEA_sens[,2,,])
mCEA_sens2 <- melt(CEA_sens[,3,,]) #new version comp with simmilar coverage
mCEA_sens2b <- dcast(mCEA_sens2, formula = Var3+Var1~ Var2 )
mCEA_sens2b$incremental_costs <- mCEA_sens2b$incremental_costs / 1e6 #costs in mio and qaly in tousands
for(i in 1:length(vaccstrat)){
  mCEA_sens2b$cut_cost_50[ mCEA_sens2b$Var1==vaccstrat[i] ] <- cut( mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],3 ],  quantile(mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],3 ], probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE) )
  mCEA_sens2b$cut_cost_95[ mCEA_sens2b$Var1==vaccstrat[i] ] <- cut( mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],3 ],  quantile(mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],3 ],  probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE) )
  mCEA_sens2b$cut_qaly_50[ mCEA_sens2b$Var1==vaccstrat[i] ] <- cut( mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],4 ], quantile(mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],4 ],probs = c(0.25,0.75) , na.rm = TRUE, names=TRUE) ) 
  mCEA_sens2b$cut_qaly_95[ mCEA_sens2b$Var1==vaccstrat[i] ] <- cut( mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],4 ], quantile(mCEA_sens2b[ mCEA_sens2b$Var1==vaccstrat[i],4 ],probs = c(0.025,0.975) , na.rm = TRUE, names=TRUE) ) 
}
# mCEA_sens2b$cut_costandqaly_quant <- mCEA_sens2b #temp if only one sensitivity scenario
# mCEA_sens2b$cut_costandqaly_quant <- ifelse(!is.na(mCEA_sens2b$cut_cost_50) & !is.na(mCEA_sens2b$cut_qaly_50), "50%", 
#                                             ifelse( !is.na(mCEA_sens2b$cut_cost_95) & !is.na(mCEA_sens2b$cut_qaly_95), "95%", NA) )
mCEA_sens2b$cut_costandqaly_IQR <- ifelse(!is.na(mCEA_sens2b$cut_cost_50) & !is.na(mCEA_sens2b$cut_qaly_50), "50%", NA) 


summary(mCEA_sens1b)
summary(mCEA_sens2b)
# mCEA_sens1b[mCEA_sens1b$Var1=="V9v_56",]
mCEA_sens2b[mCEA_sens2b$Var1=="V4v_56",]

col3 <- colorRampPalette(c("darkblue", "lightblue"))(length(diffv))
col4 <- colorRampPalette(c("darkgreen", "lightgreen"))(length(diffv))
col1 <- c( col3, col4)

#### Compared to baseline
mCEA_sens1b$ICERs <- (mCEA_sens1b$incremental_costs/ mCEA_sens1b$incremental_QALY) 

summary(mCEA_sens1b)
summary( mCEA_sens1b[which(mCEA_sens1b$ICERs> (79104/1e6) ),] )
summary( mCEA_sens1b[which(mCEA_sens1b$ICERs< 0 ),] )

summary( mCEA_sens1b[which(mCEA_sens1b$ICERs> (79104) ),] )
summary(mCEA_sens1b$ICERs[mCEA_sens1b$Var1=="V4v_80"]) *1e6
summary(mCEA_sens1b$ICERs[mCEA_sens1b$Var1=="V9v_56"]) *1e6
summary(mCEA_sens1b$ICERs[mCEA_sens1b$Var1=="V9v_80"]) *1e6


# mCEA_sens1b_sub95 <- mCEA_sens1b[!is.na(mCEA_sens1b$cut_costandqaly_quant),]
mCEA_sens1b_sub50 <- mCEA_sens1b[!is.na(mCEA_sens1b$cut_costandqaly_IQR),]

mCEA_sens1b_sub50 <- mCEA_sens1b[!is.na(mCEA_sens1b$cut_costandqaly_IQR),]

# mCEA_sens1b_sub50 <- mCEA_sens1b_sub50[mCEA_sens1b_sub50$Var1=="V9v_56",]

require(plyr)
require(ggplot2)
# Convex hulls.
df <- mCEA_sens1b_sub50
find_hull <- function(mCEA_sens1b_sub50) mCEA_sens1b_sub50[chull(mCEA_sens1b_sub50$incremental_QALY, mCEA_sens1b_sub50$incremental_costs), ]
hulls <- ddply(df, "Var1", find_hull)

pdf(paste(plotdir, "Cost-effectiveness_plane_COMPvac", Sys.Date(), ".pdf", sep=""), width=6, height=4)

ggplot(mCEA_sens1b_sub50, aes(x = incremental_QALY, y = incremental_costs, fill =  Var1 , colour=Var1 )) + 
  geom_point(  ) +  
  geom_polygon(data=hulls, alpha=.2)+
  scale_color_manual(name="Strategies" , values= col1[3:4])+
  scale_fill_manual(name="Strategies" , values= col1[3:4])+
  xlab("Incremental QALYs") + ylab("Incremental cost (in mio)")+
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  theme_bw()+
  scale_x_continuous( limits = c(-max(mCEA_sens1b_sub50$incremental_QALY)/3, max(mCEA_sens1b_sub50$incremental_QALY)+10 ) )  +
  scale_y_continuous(limits = c(min(mCEA_sens1b_sub50$incremental_costs)-1, max(mCEA_sens1b_sub50$incremental_costs)+50 ) ) +
  geom_abline(slope = 79104/1e6, linetype = "dashed")+ #79104 is the GDB per inhab in 2017 (OFS 2018)
  ggtitle("Nonavalent compared with quadrivalent vaccination")+
  theme(plot.title = element_text(size = 14, face = "bold"))

dev.off()

#### Compared to no vacc
mCEA_sens2b$ICERs <- (mCEA_sens2b$incremental_costs/ mCEA_sens2b$incremental_QALY) 
summary(mCEA_sens2b)
summary( mCEA_sens2b[which(mCEA_sens2b$ICERs> (79104/1e6) ),] )

summary( mCEA_sens2b[which(mCEA_sens2b$ICERs> (79104/1e6) ),] )
summary( mCEA_sens2b[which(mCEA_sens2b$ICERs< 0 ),] )

summary( mCEA_sens2b[which(mCEA_sens2b$ICERs> (79104/1e6) ),] )
summary(mCEA_sens2b$ICERs[mCEA_sens2b$Var1=="V4v_56"]) *1e6
summary(mCEA_sens2b$ICERs[mCEA_sens2b$Var1=="V4v_80"]) *1e6
summary(mCEA_sens2b$ICERs[mCEA_sens2b$Var1=="V9v_56"]) *1e6
summary(mCEA_sens2b$ICERs[mCEA_sens2b$Var1=="V9v_80"]) *1e6

mCEA_sens2b_sub95 <- mCEA_sens2b
# mCEA_sens2b_sub95 <- mCEA_sens2b[!is.na(mCEA_sens2b$cut_costandqaly_quant),]
mCEA_sens2b_sub50 <- mCEA_sens2b[!is.na(mCEA_sens2b$cut_costandqaly_IQR),]

# Convex hulls.
df <- mCEA_sens2b_sub50
find_hull <- function(mCEA_sens2b_sub50) mCEA_sens2b_sub50[chull(mCEA_sens2b_sub50$incremental_QALY, mCEA_sens2b_sub50$incremental_costs), ]
hulls <- ddply(df, "Var1", find_hull)

pdf(paste(plotdir, "Cost-effectiveness_plane_COMPnov", Sys.Date(), ".pdf", sep=""), width=6, height=4)

ggplot(mCEA_sens2b_sub50, aes(x = incremental_QALY, y = incremental_costs, col =  Var1, fill=Var1  )) + 
  geom_point( ) +
  geom_polygon(data=hulls, alpha=.2)+
  scale_color_manual(name="Strategies" , values= col1)+
  scale_fill_manual(name="Strategies" , values= col1[1:4])+
  xlab("Incremental QALYs") + ylab("Incremental cost (in mio)")+
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  theme_bw()+
  scale_x_continuous( limits = c(-max(mCEA_sens2b_sub50$incremental_QALY)/3, max(mCEA_sens2b_sub50$incremental_QALY)+10 ) )  +
  scale_y_continuous(limits = c(min(mCEA_sens2b_sub50$incremental_costs)-100, max(mCEA_sens2b_sub50$incremental_costs)+50 ) ) +
  geom_abline(slope = 79104/1e6, linetype = "dashed") +#79104 is the GDB per inhab in 2017 (OFS 2018)
  ggtitle("All strategies compared with no vaccination") +
  theme(plot.title = element_text(size = 14, face = "bold"))

dev.off()

unique(mCEA_sens2b_sub50$Var3)

############## Univariate analysis ####################################

############### ICER univariate sensitivity analysis 


# 1st compared to baseline vacc
#function to get ICER
#use mean estimates of the paramsets_CQok

disease_cost <- 1
utilities <- 1

#matrix for compared to baseline and one for compared to no vacc
namcol <- c("min", "base", "max")
namrow <- c("discounting", "treatement costs",  "vaccine costs", "vaccine price diff.", "utilities weights")

ICER_sens <- array(0, dim=c( 3, 3,  length(vaccstrat), length(namcol), length(namrow) ), 
                   dimnames = list( c("med", "L_IQR", "U_IQR") , c("comp56","comp80", "no vacc."),vaccstrat , namcol, namrow ))

ICER_fix_2 <- array(0, dim=c(lsens, length(vaccstrat), 3), 
                    dimnames = list(1:lsens, vaccstrat ,  c("ICER_compared_56","ICER_compared_80" ,"ICER_compared_to_novacc"  ) ))

#a discount rateps
difv <- c(0.0,0.03,0.06)

minmeanmax <- c(1,4,6)

pb1 <- txtProgressBar(min = 0, max = length(minmeanmax), style = 3)
pb2 <- txtProgressBar(min = 0, max = length(1:lsens), style = 3)

for(ii in 1:length(minmeanmax) ){
  source("01_params_script_costsqaly_means.R")
  disc <- difv[ii]
  # disc <- summary(paramsets_CQok[,which(ptv_cq=="disc")])[minmeanmax[ii]]
  for(nbsens in 1:lsens){
    ICER_fix_2[nbsens,,] <- icerf(dfinal3[nbsens,,,],nbsens,C4v, C9v, disc, disease_cost, utilities)
    setTxtProgressBar(pb2, nbsens)
  }
  for(k in 1:length(vaccstrat)){
    ICER_sens[1,1,k,ii,1] <-  summary(ICER_fix_2[,k,1])[3]
    ICER_sens[2,1,k,ii,1] <-  summary(ICER_fix_2[,k,1])[2]
    ICER_sens[3,1,k,ii,1] <-  summary(ICER_fix_2[,k,1])[5]
    
    ICER_sens[1,2,k,ii,1] <-  summary(ICER_fix_2[,k,2])[3]
    ICER_sens[2,2,k,ii,1] <-  summary(ICER_fix_2[,k,2])[2]
    ICER_sens[3,2,k,ii,1] <-  summary(ICER_fix_2[,k,2])[5]
    
    ICER_sens[1,3,k,ii,1] <-  summary(ICER_fix_2[,k,3])[3]
    ICER_sens[2,3,k,ii,1] <-  summary(ICER_fix_2[,k,3])[2]
    ICER_sens[3,3,k,ii,1] <-  summary(ICER_fix_2[,k,3])[5]
  }
  setTxtProgressBar(pb1, ii) }

ICER_sens[1,1,,2,]

ICER_sens[1,1,1,1,1]
ICER_sens[,1,2,,2]

#b treatement costs
IQRmean <- c(2,4,5) #takes the lower IQR, mean and higher IQR of the lognormal distribution
for(ii in 1:length(IQRmean) ){
  source("01_params_script_costsqaly_means.R")
  Ct[1,1] <- summary(paramsets_CQok[,which(ptv_cq=="Ct_1")])[IQRmean[ii]]
  Ct[1,2] <- summary(paramsets_CQok[,which(ptv_cq=="Ct_2")])[IQRmean[ii]]
  Ct[1,3] <- summary(paramsets_CQok[,which(ptv_cq=="Ct_3")])[IQRmean[ii]]
  Ct[1,4] <- summary(paramsets_CQok[,which(ptv_cq=="Ct_4")])[IQRmean[ii]]
  Ct[1,5] <- summary(paramsets_CQok[,which(ptv_cq=="Ct_5")])[IQRmean[ii]]
  
  Cpt <- summary(paramsets_CQok[,which(ptv_cq=="Cpt")])[IQRmean[ii]]
  Cfp <- summary(paramsets_CQok[,which(ptv_cq=="Cfp")])[IQRmean[ii]]
  # Cme <- summary(paramsets_CQok[,which(ptv_cq=="Cme")])[IQRmean[ii]]
  # Cse <- summary(paramsets_CQok[,which(ptv_cq=="Cse")])[IQRmean[ii]]
  for(nbsens in 1:lsens){
    ICER_fix_2[nbsens,,] <- icerf(dfinal3[nbsens,,,],nbsens,C4v, C9v, disc, disease_cost, utilities)
    setTxtProgressBar(pb2, nbsens)
  }
  for(k in 1:length(vaccstrat)){
    ICER_sens[1,1,k,ii,2] <-  summary(ICER_fix_2[,k,1])[3]
    ICER_sens[2,1,k,ii,2] <-  summary(ICER_fix_2[,k,1])[2]
    ICER_sens[3,1,k,ii,2] <-  summary(ICER_fix_2[,k,1])[5]
    
    ICER_sens[1,2,k,ii,2] <-  summary(ICER_fix_2[,k,2])[3]
    ICER_sens[2,2,k,ii,2] <-  summary(ICER_fix_2[,k,2])[2]
    ICER_sens[3,2,k,ii,2] <-  summary(ICER_fix_2[,k,2])[5]
    
    ICER_sens[1,3,k,ii,2] <-  summary(ICER_fix_2[,k,3])[3]
    ICER_sens[2,3,k,ii,2] <-  summary(ICER_fix_2[,k,3])[2]
    ICER_sens[3,3,k,ii,2] <-  summary(ICER_fix_2[,k,3])[5]
  }
  setTxtProgressBar(pb1, i) }


#c vaccination price 
for(ii in 1:length(IQRmean) ){
  source("01_params_script_costsqaly_means.R")
  pricediff <- mean(paramsets_CQok[,which(ptv_cq=="C9v")] - paramsets_CQok[,which(ptv_cq=="C4v")])
  C4v <-  summary(paramsets_CQok[,which(ptv_cq=="C4v")])[IQRmean[ii]]
  C9v <- C4v + pricediff
  for(nbsens in 1:lsens){
    ICER_fix_2[nbsens,,] <- icerf(dfinal3[nbsens,,,],nbsens,C4v, C9v, disc, disease_cost, utilities)
    setTxtProgressBar(pb2, nbsens)
  }
  for(k in 1:length(vaccstrat)){
    ICER_sens[1,1,k,ii,3] <-  summary(ICER_fix_2[,k,1])[3]
    ICER_sens[2,1,k,ii,3] <-  summary(ICER_fix_2[,k,1])[2]
    ICER_sens[3,1,k,ii,3] <-  summary(ICER_fix_2[,k,1])[5]
    
    ICER_sens[1,2,k,ii,3] <-  summary(ICER_fix_2[,k,2])[3]
    ICER_sens[2,2,k,ii,3] <-  summary(ICER_fix_2[,k,2])[2]
    ICER_sens[3,2,k,ii,3] <-  summary(ICER_fix_2[,k,2])[5]
    
    ICER_sens[1,3,k,ii,3] <-  summary(ICER_fix_2[,k,3])[3]
    ICER_sens[2,3,k,ii,3] <-  summary(ICER_fix_2[,k,3])[2]
    ICER_sens[3,3,k,ii,3] <-  summary(ICER_fix_2[,k,3])[5]
  }
  setTxtProgressBar(pb1, ii) }


#d vaccination price difference
for(ii in 1:length(IQRmean) ){
  source("01_params_script_costsqaly_means.R")
  pricediff <- runif( lsens,  vpd*0.1, vpd*10 )
  C9v <- C4v + summary(pricediff)[IQRmean[ii]]
  for(nbsens in 1:lsens){
    ICER_fix_2[nbsens,,] <- icerf(dfinal3[nbsens,,,],nbsens,C4v, C9v, disc, disease_cost, utilities)
    setTxtProgressBar(pb2, nbsens)
  }
  for(k in 1:length(vaccstrat)){
    ICER_sens[1,1,k,ii,4] <-  summary(ICER_fix_2[,k,1])[3]
    ICER_sens[2,1,k,ii,4] <-  summary(ICER_fix_2[,k,1])[2]
    ICER_sens[3,1,k,ii,4] <-  summary(ICER_fix_2[,k,1])[5]
    
    ICER_sens[1,2,k,ii,4] <-  summary(ICER_fix_2[,k,2])[3]
    ICER_sens[2,2,k,ii,4] <-  summary(ICER_fix_2[,k,2])[2]
    ICER_sens[3,2,k,ii,4] <-  summary(ICER_fix_2[,k,2])[5]
    
    ICER_sens[1,3,k,ii,4] <-  summary(ICER_fix_2[,k,3])[3]
    ICER_sens[2,3,k,ii,4] <-  summary(ICER_fix_2[,k,3])[2]
    ICER_sens[3,3,k,ii,4] <-  summary(ICER_fix_2[,k,3])[5]
  }
  setTxtProgressBar(pb1, i) }

##### 
#c utilies weights
for(ii in 1:length(minmeanmax) ){
  source("01_params_script_costsqaly_means.R")
  q[1,1] <-  summary(paramsets_CQok[,which(ptv_cq=="q_1")])[minmeanmax[ii]]
  q[1,2] <-  summary(paramsets_CQok[,which(ptv_cq=="q_2")])[minmeanmax[ii]]
  q[1,3] <-  summary(paramsets_CQok[,which(ptv_cq=="q_3")])[minmeanmax[ii]]
  q[1,4] <-  summary(paramsets_CQok[,which(ptv_cq=="q_4")])[minmeanmax[ii]]
  q[1,5] <-  summary(paramsets_CQok[,which(ptv_cq=="q_5")])[minmeanmax[ii]]
  # q <- ifelse( (q * utilities) > 1, 1, q * utilities)
  for(nbsens in 1:lsens){
    ICER_fix_2[nbsens,,] <- icerf(dfinal3[nbsens,,,],nbsens,C4v, C9v, disc, disease_cost, utilities)
    setTxtProgressBar(pb2, nbsens)
  }
  for(k in 1:length(vaccstrat)){
    ICER_sens[1,1,k,ii,5] <-  summary(ICER_fix_2[,k,1])[3]
    ICER_sens[2,1,k,ii,5] <-  summary(ICER_fix_2[,k,1])[2]
    ICER_sens[3,1,k,ii,5] <-  summary(ICER_fix_2[,k,1])[5]
    
    ICER_sens[1,2,k,ii,5] <-  summary(ICER_fix_2[,k,2])[3]
    ICER_sens[2,2,k,ii,5] <-  summary(ICER_fix_2[,k,2])[2]
    ICER_sens[3,2,k,ii,5] <-  summary(ICER_fix_2[,k,2])[5]
    
    ICER_sens[1,3,k,ii,5] <-  summary(ICER_fix_2[,k,3])[3]
    ICER_sens[2,3,k,ii,5] <-  summary(ICER_fix_2[,k,3])[2]
    ICER_sens[3,3,k,ii,5] <-  summary(ICER_fix_2[,k,3])[5]
  }
  setTxtProgressBar(pb1, ii) }

date_save <- gsub("-", "_", Sys.Date()) 

save(ICER_sens, file=paste(  "ICER_sens" , date_save, ".RData", sep="") )

### PLOT IT
require("plotrix")

# compared to baseline

# namrow1 <- namrow[-5]
namrow1 <- namrow

# nbparam <- (length(namrow)-1)
nbparam <- (length(namrow))

sepg <- 1.5
yax <- c(1:nbparam, (nbparam+sepg) : ( 2*nbparam+sepg-1 )  )
# yax <- c(1:nbparam, (nbparam+sepg) : ( 2*nbparam+sepg-1 )  ,
#          ( 2*nbparam+2*sepg-1 ) : ( 3*nbparam+2*sepg-2 )  )

yax <- rev( yax )  -0.5
sepg1 <-  c(1,(nbparam+sepg)   ) -0.5
# sepg1 <-  c(1,(nbparam+sepg) ,( 2*nbparam+2*sepg-1 ) ,( 3*nbparam+3*sepg-2 )  ) -0.5

pdf(paste(plotdir, "ICER_univ_sens_v2", Sys.Date(), ".pdf", sep=""), width=6, height=6)
par(mar=c(4,7,2,10))
plot(NA, NA, pch=19, frame=FALSE, main="Nonavalent compared with quadrivalent vaccination",
     ylab="", xlab="ICER", yaxt="n", ylim=c(0,max(yax)+1.5),
     xlim=c(min(ICER_sens[2,1,3,,-5], na.rm=TRUE),max(ICER_sens[3,1,3,,-5]+1700, na.rm=TRUE) ) )

for(j in 1:nbparam){
  for(i in 1:2){
    points(ICER_sens[1,i,2+i,2,j ], yax[(i*nbparam)-(nbparam-j )], pch=19+j, bg="white",cex=0.75)
    segments(x0 =  ICER_sens[1,i,2+i,1,j ],   y0= yax[(i*nbparam)-(nbparam-j )],
             x1 = ICER_sens[1,i,2+i,3,j ],    y1= yax[(i*nbparam)-(nbparam-j )]) } }

staxlab(4,at=yax,labels=namrow1, srt=0, cex=0.75, font=1, adj=-0)
abline(h = sepg1-0.75, col="grey", lty=3, lwd=0.7)
abline(v = 0, col="black", lty=1, lwd=0.2)
staxlab(2,at=(sepg1+(nbparam/2))[-4] -2.5,labels=c("80% V9v", "56% V9v"), srt=270, cex=1, font=2)
# staxlab(4,at=(sepg1+(nbparam/2))[-4] -1.8,labels=c("80% V9v", "56% V9v","80% V4v"), srt=270, cex=1, font=2)
dev.off()

# compared to no vacc
#nbparam <- (length(namrow)-1)
# namrow1 <- namrow[-5]
namrow1 <- namrow

# nbparam <- (length(namrow)-1)
nbparam <- (length(namrow))

sepg <- 1.5
yax2 <- c(1:nbparam, (nbparam+sepg) : ( 2*nbparam+sepg-1 )  ,
          ( 2*nbparam+2*sepg-1 ) : ( 3*nbparam+2*sepg-2 ),
          ( 3*nbparam+3*sepg-2 ):  ( 4*nbparam+3*sepg-3 ) )  
yax2 <- rev( yax2 )  -0.5

sepg1 <-  c(1,(nbparam+sepg) ,( 2*nbparam+2*sepg-1 ) ,( 3*nbparam+3*sepg-2 ), ( 4*nbparam+4*sepg-3 )  ) -0.5

pdf(paste(plotdir, "ICER_univ_sens_comp_nov_v2_with_UT", Sys.Date(), ".pdf", sep=""), width=6, height=6)
par(mar=c(4,7,2,10))

plot(NA, NA, pch=19, frame=FALSE, main="All strategies compared with no vaccination",
     ylab="", xlab="ICER", yaxt="n", ylim=c(0,max(yax2)+1.5),
     xlim=c(min(ICER_sens[2,3,,,-5], na.rm=TRUE),max(ICER_sens[3,3,,,-5]+1, na.rm=TRUE) ) )
for(j in 1:nbparam){
  for(i in 1:4){
    points(ICER_sens[1,3,i,2,j ], yax2[(i*nbparam)-(nbparam-j )], pch=19+j, bg="white",cex=0.75)
    segments(x0 =  ICER_sens[1,3,i,1,j ],   y0= yax2[(i*nbparam)-(nbparam-j )],
             x1 = ICER_sens[1,3,i,3,j ],    y1= yax2[(i*nbparam)-(nbparam-j )]) } }
# for(j in 1:nbparam){
#   for(i in 1:4){
#     points(ICER_sens[1,2,i,2,j ], yax2[(i*nbparam)-(nbparam-j )], pch=19+j, bg="white",cex=0.75)
#     segments(x0 =  ICER_sens[1,2,i,1,j ],   y0= yax2[(i*nbparam)-(nbparam-j )],
#              x1 = ICER_sens[1,2,i,3,j ],    y1= yax2[(i*nbparam)-(nbparam-j )]) } }

staxlab(4,at=yax2,labels=namrow1, srt=0, cex=0.75, font=1, adj= -0)
abline(h = sepg1-0.75, col="grey", lty=3, lwd=0.7)
abline(v = 0, col="black", lty=1, lwd=0.2)
staxlab(2,at=(sepg1+(nbparam/2))[-5] -2.6,labels=c("80% V9v", "56% V9v","80% V4v", "56% V4v"), srt=270, cex=1, font=2)

dev.off()








  
  
  