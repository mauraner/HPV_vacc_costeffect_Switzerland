# 01_params_script_baselineparam
# MR, May 2019

# Retrieves baseline parameters for the main HPV model transmission dynamics

#--------------------------- directory to take data from
dirl <- "tab_params/"

#--------------------------- Model dimensions ------------------------------

sex <- c("f", "m")

riskgroups <- c("l", "h")

agegroups <- c("10-14", "15-19", "20-24","25-29" ,"30-39", "40-64", "65-99") #max 72 for natsal3, and 99 for shs

comp <- c("S", "I", "CIN2", "CIN3", 
          "ULCC", "DLCC", "URCC", "DRCC",  "UDCC", "DDCC", "DCC", "R", "V")

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

dis_outcome <- c("cervix", "anus", "penile", "oropharynx", "warts")

vaccines <- c("no","quadrivalent", "nonavalent")



#--------------------------- FIXED params ------------------------------

# Vaccination  **********
#vaccination parameters are set to 0, this can be changed in the correspondinng script before running the model, according to prevention strategy choosen

# proportion young women vaccinated, has to be defined per type: w==1 | w==2 | w== 11 | w==12 for quadrivalent,
#w==1 | w==2 |w==3 |w==4 | w==7 |w==9 |w==10 |  w== 11 | w==12 for nonavalent, non-vaccine types are : 5,6,8,13

#vaccination coverage rates p[sex,riskgrous,agegroups] for quadrivalent vaccine
pv <- array(0, dim=c(length(sex),length(hpvtype), length(vaccines)), dimnames = list(sex,hpvtype, vaccines )  )

vs <- 1 # vaccine strategy 1 is no vaccination, 2 is quadrivalent 3 is nonavalent vaccination

# population **********

#population size N[sex,riskgroup,agegroup], the sum of all the population (including deaths) is 1
N <- array(1, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups,
                                                                                       agegroups)  )
#sex ratio
sr <- 0.5

N[1,,] <- N[1,,]*sr 
N[2,,] <- N[2,,]*sr 

#risk groups % in high or low sexual activity groups. 
# This is retrieved in the scripts: 01_params_compl_partnerchangerate1_CHsexbehav and 01_params_compl_partnerchangerate2_MLEactivitygroups
# the resulting data is saved as sbehav_ch3.csv
sexbehavdata <- paste(dirl, "sbehav_ch3.csv", sep="")
sbehav <- as.data.frame(read.csv(sexbehavdata) )

rgr <- c(sbehav[1,2],1-sbehav[1,2]) # % in low or high sexual activity group

N[,1,] <- N[,1,]*rgr[1]
N[,2,] <- N[,2,]*rgr[2]


#makes age groups and adapt size according to number of year spent in each group 
agr1 <- c(14-10 ,(19-15) , 24-20,  29-25, 39-30, 64-40,99-65)
agr1 <- agr1 +1
agr <- agr1/ mean(agr1)

agr <- agr* (1/ length(agegroups) )


for( i in 1:length(agegroups) ){
  N[,,i] <- N[,,i]*agr[i]
}

#rate of entering/exiting system, according to the min and max age included
mu <- 1/ (100-10)

#rates of changes between agegroups
g <- mu/agr

g <- matrix(g, nrow=1, ncol=length(agegroups))
colnames(g) <- agegroups

#number of contacts C[sex, riskgroups, age]
C <- array(0, dim=c(length(sex),length(riskgroups),length(agegroups)), dimnames = list(sex,riskgroups, agegroups)  )

# sexbehavdata <- "sbehav_ch3.csv"
sbehav <- as.data.frame(read.csv(sexbehavdata) )
rownames(sbehav) <- sbehav[,1]
sbehav1 <- sbehav[-1,2:3]

for( i in 1: (length(agegroups)-1) ){
  C[,1,i+1] <- sbehav1[c( (i*2) -1 , (i*2) ),1]
  C[,2,i+1] <- sbehav1[c( ( ( (i*2) -1) + (length(agegroups)-1)*2), ( (i*2) + (length(agegroups)-1)*2) ), 1]
}

C[,,1] <- 1e-31


#sexual mixing between risk groups (eps_r) and age groups (eps_a), assortativity (0= random, 1=assortative), 
# eps_r is based on the litterature, and eps_a was retrieved from information hold in Natsal'3 dataset (01_params_compl_agegroupassortativity.R)
eps_r <- 0.5
eps_a <- 0.5242733

# m is the proportion who can change sexual activity group over time, redistributed accoring to the proportion in each group
# value 1 means that people can change, however, since low activity group is the biggest, they mainly stay in that group
# based on litterature
m <- 1


#--------------------------- PER HPV TYPE VARYING PARAMETERS ------------------------------

hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "otherHR")

# ********** infection dynamics **********

#transmissibility probablity
betap <- array(0, dim=c(length(sex),length(hpvtype) ), dimnames = list(sex,hpvtype) )

#scaled to fit to prevalence data. Here is only an example from previous optimisation. Has to be called before running the model
# can be adapted to women to men and men to women, however, we assume the same for both sexes
# betap[1,1:13] <-c(0.95, 0.93, 0.62, 0.89, 0.66, 0.61, 0.68, NA, 0.63, 0.65,NA, NA,  0.89)
betap[1,1:13] <- c(0.9885419, 0.8276534, 0.8355821, 0.7789747 ,0.7501864, 0.7714783 ,0.7890471, 
                   0.8007531, 0.8174955, 0.7749328, 0.549997,0.5499666 ,0.8827871)
betap[2,1:13] <- betap[1,1:13] 

#Proportion of infections that develop HPV-type specific antibodies #from Durham and for 6 11 and other it is the mean value
alpha <- array(0, dim=c(length(hpvtype) ), dimnames = list(hpvtype) )

temphpv <- binom.test((7+1+68+16+47+11+8+19+11+10+8+11+17+9+9),(11+1+120+39+83+27+9+40+29+65+22+19+41+20+22))  #for hpv51, 11 and OHR, all HPV types from the paper were pooled together to get an average
temphpv2 <- as.numeric(temphpv$estimate)*100

alpha[1:13] <- c(56.7, 41.0, 56.6, 40.7, 88.9, 47.5, 37.9, temphpv2, 15.4, 57.9, 63.6, temphpv2, temphpv2 ) #for hpv51, 11 and OHR, all HPV types from the paper were pooled together to get an average
alpha[1:13] <- alpha[1:13]/100


#--------------------------- PER DISEASE OUTCOME VARYING ------------------------------
# It was meant to include different HPV disease outcomes. However, in the manuscript, only cervix is used

dis_outcome <- c("cervix", "anus", "penile", "oropharynx", "warts")

#rate of progression to different stage: pri1 = from infectious status to CIN2, pri2 = from infectious to CIN3
pri1 <- array(0, dim=c(length(sex),length(hpvtype), length(dis_outcome) ), dimnames = list(sex,hpvtype, dis_outcome) )
pri2 <- array(0, dim=c(length(sex),length(hpvtype), length(dis_outcome) ), dimnames = list(sex,hpvtype, dis_outcome) )


#from Jaisamarm 2013, take the HR from the multivariate analysis compared to the probability in the non-oncogenic types
pri1[1,,1] <- -log(1- ( (42/4824) * c(9.25,3.56,5.09,9.14, rep(2.63,2), 3.64, rep(2.63,3),0,0, 2.63) ))/2

pri2[1,,1] <- -log(1- ( (7/4824) * c(20.93,4.74,7.82,20.47, rep(3.51,2), 4.45, rep(3.51,3),0,0, 3.51) ))/2


# progression from CIN2 to CIN3, From Durham reference: Moscicki A-B et al. 2010, proportion by year 3 (15%, 95%CI:9%-26%)
pr23 <- array(0.05417298, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) )
pr23[c(11,12),1:4] <- 0 # for low risk types for other outcomes than genital warts
pr23[-c(11,12),5] <- 0

# progression from CIN3 to cancer is scaled in incidence data. This is just an example based on previous calibrations, data has to be loaded before running the model
# pr3c <- c(0.040, 0.080 ,0.095, 0.005, 0.020 ,0.175,0.105, 0.0, 0.045, 0.040, 0,0, 0.010)
pr3c <- c(0.036116553, 0.042435846, 0.009402965, 0.016722004, 0.057802549, 0.016278879,
          0.063148826, 0.017084773, 0.020716105,0.021775352, 0.000000000, 0.000000000, 0.050164198)
names(pr3c) <- hpvtype
# progression from CIN2 to cancer is dependent on progression CIN3 to cancer, which is scaled to cancer incidence (pr2c = 0.2*pr3c), From Campos et al. 2014
# This is just an example, data has to be loaded before running the model
# pr2c <- c(0.040, 0.080 ,0.095, 0.005, 0.020 ,0.175,0.105, 0.0, 0.045, 0.040, 0,0, 0.010)*0.2
pr2c <- 0.2* pr3c
names(pr2c) <- hpvtype


#rate of regression from CIN3 and CIN2 to susceptibles # from Moscicki 2010 for CIN2 and McCredie for CIN3
pr3s <- array(0, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) )
pr2s <- array(0, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) ) 

pr3s <- array(0.1, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) )
pr3s[c(11,12),2:5] <- 0

pr2s <- array(0.5047092, dim=c(length(hpvtype), length(dis_outcome) ), dimnames = list(hpvtype, dis_outcome) ) #Moscicki 2010
pr2s[1:2,1] <- 0.2669108

#This was assumed to be an unnecessary parameters since we scale beta and pr2c/pr3c. Also, the estimates were not very relieable in my opinon (McCredie et al 2008)
pr32 <- 0

#progression localised to regional and regional to distant, from Campos et al (2014), slightly different from Durham's
prlcrc <- 0.2424 #
prrcdc <- 0.3038

#relative risk of re-infection following clearance, from Castellsagué et al. 2014
tau <- array(0, dim=c(length(hpvtype), length(dis_outcome) ), 
             dimnames = list(  hpvtype, dis_outcome) )

tau[,1:5] <- c(0.64, rep(0.96,12) ) #only known for hpv16 and 18. All other types are considered to be like hpv18
# CI: 0.64 (0.53-0.78), 0.94 (0.75-1.19)

#Vaccine efficacy #hpv16:0.91-0.96, 18:0.91-0.97): from Lu et al 2011. 
#for the other hpv types: 0.967 CI 0.809-0.998) from Joura et al 2015
ve <- rep(0, length(hpvtype))
ve <- c(0.94, 0.95, 0.967, 0.967, 0, 0,0.967,0,0.967, 0.967, 0.967, 0.967, 0)
names(ve) <- hpvtype

#screening rates per age group and disease (only women)
#take % of women who had pap smear in the last year from SHS 2012, retrieved from SHS 2012 through: 01_params_compl_screeningCH.R, saved under: sceeningrates_2012SHS_V2.csv
s1  <- read.csv( paste(dirl, "sceeningrates_2012SHS_v2.csv", sep="") )
s1 <- as.data.frame(s1[,2:13])

s <- array(0, dim=c(length(sex), length(agegroups), length(dis_outcome) ), 
           dimnames = list( sex, agegroups, dis_outcome) )
#screening only occurs for cervix outcomes and only in women # for the first age group there was no indication, I took half of the second age group
# s[1,,1] <-  c( s1[1,6] /2 , s1[1:5,6] )
s[1,,1] <-  c( 0 ,0, s1[c(1:5),6] )

#transform proportions in rate.
s[1,,1] <- (-log(1-s[1,,1]) )/1

# mt is a factor which accounts for the increased effetctiveness of screening in case of multitype infections
# it is added in the model as s+mt*s, mt being the proportion of multitype infections, retreived from the CIN3+plus study. With this,
# we account for the fact that the type might be screening during the "screening" of another type

mt <- 0.128

#screening sensitivity
ssens  <- 0.587 #from Bigras 2005, CI48.6-68.2

#screening specificity
sspec <- 0.969 #from Bigras 2005, CI 96.6-97.2

#--------------------------- PER DISEASE STAGE SPECIFIC (for cervical cancer) ------------------------------

#stage specific probablity of diagn.(?), from Campos (2014) Web appendix, table 1
dstages <- c("cin2", "cin3", "local", "regional", "distant")
dstages[3:5]

z <- array(0, dim=c(length(dstages[3:5]), length(dis_outcome) ), 
           dimnames = list( dstages[3:5],  dis_outcome) )
z[,1:5] <- c(0.2106379, 0.9160948,2.302646)


#stage specific probability of cure with treatement, from Durham who cites: Elbasha EH et al. 2007 # cure rate
pi <- array(0, dim=c(length(dstages[3:5]), length(dis_outcome) ), 
            dimnames = list( dstages[3:5],  dis_outcome) )
pi[,1:5] <-  c(0.92, 0.55, 0.17)

#mortality of undiagnosed and diagnosed cancers (stage specific), from Campos et al. (2014), slightly different from Durham but very close
um <- c(0.01921,0.11455,0.35685)
dm <- c(0.01080486,0.04327795,0.09154833)


#infection clearance rate (1/gamma = duration of the infectiousness), or disease clearance rate
#gamma[infection clearance, CIN2 clearance, CIN3 clearane]
gamma <- array(0, dim=c(length(sex),length(hpvtype)), dimnames = list(sex,hpvtype) ) 

# per hpv type from Jaisamararn et al (2013). Different from Durham # median duration, no multivariate analysis
# gamma[1,] <- c(0.49, 0.66, 0.57, 0.66, 0.69,0.69, 0.66, 0.69, 0.69 ,0.69, 0.69, 0.69, 0.69 ) #these are Durhams
gamma_a1 <- c(17.11, 11.84, 13.8, 12.00, 11.77,11.77, 11.48, 11.77, 11.77 ,11.77, 8.26, 8.26, 11.77 )
gamma[1,] <- 1/ ( gamma_a1 /12 )
#for men , from Moreira 2014              
# gamma[2,] <- c(1.08, 1.33, rep(1.41,11 )  )
gamma_a2 <- c(7.7, 6.2, rep(6.1,11 )  )
gamma[2,] <- 1/ ( gamma_a2 /12 )

#immunity duration (1/omega= duration for becoming susceptible again)
######### New from Johnson (2012) paper because I could not retrieve the values in the Durham reference (Syrjänen 2009)#----------------------------------CHECK OK
omega <- array(0, dim=c(length(sex),length(hpvtype)), dimnames = list(sex,hpvtype) ) 

tempo <- c(0.024, 0.017, 0.018, 0.018 , 0.016, 0.021, 0.016, 0.019, 0.027, 0.017,0.020,0.017  ) # these are all available HPV HR values
mean(tempo)
omega[1,] <- c(0.024, 0.017, 0.018, 0.018 , 0.016, 0.021, 0.016, 0.019, 0.027, 0.020,mean(tempo), mean(tempo), mean(tempo)  ) 
omega[2,] <- omega[1,]

# END of parameters for infection dynamics.

#--------------------------- Initial parameters and name vectors ------------------------------

#initial params I0[compartement, sex, riskgroups, agegroups]
I0_1<- array(0, dim=c( length(comp), length(sex), length(riskgroups), length(agegroups)),
             dimnames =list(comp,sex,riskgroups,
                            agegroups))

#set initial infection proportion in the population for the different risk groups and age group
#here I assume it is the same for both age groups and for male and female

# women & men
# Suseptibles is N-I0
I0_1[1,,1,] <-  N[,1,]- 0.01*N[,1,] 
# highrisk group
I0_1[1,,2,] <- N[,2,]- 0.1*N[,2,]


# lowrisk group
I0_1[2,,1,] <- 0.01*N[,1,] 
# highrisk group
I0_1[2,,2,] <- 0.1*N[,2,]


I0 <- I0_1
I0[,1,,1] / (sr* agr[1])


##### Output column names vector
# gives column names to output matrix
compnamef <- function(compartements){
  c_name <- rep(compartements, length(sex)*length(riskgroups)*length(agegroups))
  s_name <- rep( rep(sex, each=length(compartements)) , length(riskgroups)*length(agegroups))
  r_name <- rep( rep(riskgroups, each=length(compartements)*length(sex) ), length(agegroups))
  a_name <- rep(agegroups, each=length(compartements)*length(sex)*length(riskgroups) )
  col_nam <- c("time", paste(c_name, s_name, r_name, a_name, sep="_") )
  return(col_nam)
}

col_nam <- compnamef(comp) 




#--------------------------- Age group weights ------------------------------


#age group weights according to the European Standard population (ESP 2013)
#agw <- c( sum(0.01,0.04,0.055,0.055),0.055,0.06,0.06, sum(0.065,0.07), sum(0.07,0.07,0.07,0.065,0.06), sum(0.055,0.05,0.04,0.025), 0.025 )
#names(agw) <- c("0:14",agegroups, "+85")

# from http://ec.europa.eu/eurostat/documents/3859598/5926869/KS-RA-13-028-EN.PDF/e713fa79-1add-44e8-b23d-5e8fa09b3f8f p.31, projections 2011-2020
nav <- c(1078.641, 4373.749, 5410.346, 5252.859,  5410.049, 6066.914,  6711.973,  7023.97, 7135.495, 
         7126.248, 7087.804, 6938.434, 6595.514, 6095.677, 5307.002, 4328.78, 3419.627, 2492.941,1452.548,555.307,136.12 )
names(nav) <- c("0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
                "55-59", "60-64","65-69", "70-74","75-79", "80-84" , "85-89","90-94", "95+")

nav /sum(nav)

#from NICER reference of ESP (old ESP version): http://www.nicer.org/NicerReportFiles2018/FR/methods_file/methods.htm
nav2 <- c(8, rep(7,10), 6, 5,4,3,2,1,1)
nav2 <- nav2*1000
names(nav2) <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
                 "55-59", "60-64","65-69", "70-74","75-79", "80-84" , "85+")

sum(nav2)
nav2 <- nav2/sum(nav2)

#make a vector with age group weights which corresponds to our age groups and the younger population (not included in the model)
agw <- c( sum(nav[1:3]),  nav[4:7],  sum(nav[8:9]), sum(nav[10:14]), sum(nav[15:21]) )
names(agw) <- c("0:9",agegroups)

agw <- c( sum(nav2[1:2]),  nav2[3:6],  sum(nav2[7:8]), sum(nav2[9:13]), sum(nav2[14:18]) )
names(agw) <- c("0:9",agegroups)

#Swiss age weights based on age groups used, from the Bundesamt für Statistik 2017 - data from 2016
sagw_r <- as.data.frame( read.csv( paste(dirl, "OFS_2017_pop_age.csv", sep="") ) )

totpop <- sum(sagw_r)

sagw <- c( sum(sagw_r[1:10]), sum(sagw_r[11:15]),  sum(sagw_r[16:20]), sum(sagw_r[21:25]), sum(sagw_r[26:30]),
           sum(sagw_r[31:40]),  sum(sagw_r[41:65]) ,sum(sagw_r[66:101]) ) / totpop
names(sagw) <- c("0:9",agegroups)
sum(sagw)


#--------------------------- Cervical Cancer incidence ------------------------------
# Retrieve CC incidence from NICER crude estimates from 2004-2014
# uses 01_params_compl_ccIncidenceandtyperatio.R script
totinc_cc_tot <- 6.427273 #totinc_cc

#incidence for modelled age groups (excludes <10 y.o)
totinc_cc_1099 <- totinc_cc_tot / sum(sagw[2:8])

# *** Cervical incidence per HPV type ***

# type ratio source (2010 de Sanjosé, Lance Oncol.)
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")
hpvtypprop <- hpvtype
names(hpvtypprop) <- hpvtype

notmodelledhpvtypesHR <- c(3, 5, 1,1, 44,32,15,2,3,13,4,1,1,10)
names(notmodelledhpvtypesHR) <- c("26","30", "34", "39","44","56", "59", "66","67", "68", "68b", "70", "91" , "underter")

names(hpvtypprop) <- hpvtype
hpvtyppropnb <- c(1348, 150, 69, 117, 46, 27, 80, 28, 40 , 27,  0.0,0.0, sum(notmodelledhpvtypesHR) )

hpvtypprop <- hpvtyppropnb / sum(hpvtyppropnb)

cc_inc_pertype_tot <- hpvtypprop * totinc_cc_tot # Incidence per type per 100'000 women (total pop)
cc_inc_pertype_1099 <- hpvtypprop * totinc_cc_1099 # same but for model age group (excl <10 y.o)

sum( cc_inc_pertype_tot[c(1,2,3,4,7,9,10,11,12)] ) / sum(cc_inc_pertype_tot) #prop in nonavalent vaccine
sum(hpvtyppropnb)

# write.csv(cc_inc_pertype_tot, file = "tab_params/cc_inc_pertype4e6.csv")

dirlts <- "C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/" #where to save the latex tables to


as.vector(cc_inc_pertype_tot)

# write.table(round(cc_inc_pertype_tot,2)[-c(11,12)], paste(dirlts,"cc_inc_CH_totpop.txt", sep=""), quote=FALSE, eol="&", sep=" & ", col.names = F, row.names = FALSE)


#--------------------------- HPV prevalence and type ratio for normal cytology ------------------------------

# uses 01_params_compl_ccIncidenceandtyperatio.R script, based on Guan et al 2012

#load type ratio for normal cytology
nGUAN_typerationormal <- read.csv(paste(dirl, "nGUAN_typerationormal.csv", sep=""))[,2]
#calculate type specific prevalence based on Natsal'3 estimates of HR HPV in 15-19 years olds (from 01_params_compl_hpvprevnatsal.R)
#attribute this HR prevalence to each HPV type according to the type ratio from Guan et al. and assumes that this corresponds to the 
#HPV HR prevalence in 20-24 year old Swiss women
prev_pt_2024 <- nGUAN_typerationormal * 0.2435849
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")
names(prev_pt_2024) <- hpvtype

# write.table(round(prev_pt_2024*100,2)[-c(11,12)], paste(dirlts,"hpvprev2024.txt", sep=""), quote=FALSE, eol="&", sep=" & ", col.names = F, row.names = FALSE)


