# MR, 23.01.2019

########### Parameters for the costs and qaly aspects of main HPV model


# copied from R_ce folder: CE_model_params_v7.R


#--------------------------- Parameters: ONLY COSTS ------------------------------

### COSTS
#annual discount rate
dr<- 0.03

#life-year weight (per stage: CIN2, CIN3, Local, regional, distant)
dstages <- c("cin2", "cin3", "local", "regional", "distant")
q <- array(0, dim=c(length(dstages),length(dis_outcome)), dimnames = list(dstages,dis_outcome) ) 
#cervix
q[,1] <- c(0.92,0.91,0.76,0.67,0.48)
#anus
q[,2] <- c(0.92,0.91,0.76,0.67,0.48)
#penile
q[,3] <- c(0.92,0.91,0.76,0.67,0.48)
# oropharynx
q[,4] <- c(0.92,0.91,0.76,0.67,0.48)
# warts
q[,5] <- c(0.2,0,0,0,0)


#treatment costs (per stage)
Ct <- array(0, dim=c(length(dstages),length(dis_outcome)), dimnames = list(dstages,dis_outcome) ) 
Ct[,1] <- c(1150.7, 2237.6, 20000, 20000+1890, 20000+1890+17110) # for cin2-3 and local cancer: Swiss Szucs 2008 study, then I add the difference between local and regional and betwen regional and distant


#Vaccine costs (three dose series)
# C2v <- 389
C4v <- 2* (66.6 +23.7) #Price is based on a press release from GDK. CDS from 2008
C9v <- 2*(66.6 +23.7+ (195.10-190.94 ) )  # I added the price difference between 1 V4v hpv dose and one V9v dose informed by MSD
C9v <- 2*(66.6 +23.7+ (195.10-190.94 ) )  # I added the price difference between 1 V4v hpv dose and one V9v from Durham 

vpd <- C9v-C4v

# C4v <- 2* (160 +23.7) #Price is based on a press release from GDK. CDS from 2008
# C9v <- 2*(160 +23.7+ (195.10-190.94 ) )  # I added the price difference between 1 V4v hpv dose and one V9v dose from Durham paper

#costs per negative cercival screen
Cpt <- 54.5+155.0 #from Szucs 2008 paper  # Durham: 96.42


#costs of follow-up of false positive
Cfp <- 80 #I take colposcopy/biopsy from Szucs 2008 #Durham:305.66
#cost per mild adverse effect (ae)
Cme <- 45
#QALY per mild adverse effect (ae)
Qme <-0.001
#cost per mild adverse effect (ae)
Cse <- 3280
#QALY per mild adverse effect (ae)
Qse <- 0.076

#rate of mild adverse effect
rme <- 105/100000
#rate of severe adverse effect
rse <- 9/100000




