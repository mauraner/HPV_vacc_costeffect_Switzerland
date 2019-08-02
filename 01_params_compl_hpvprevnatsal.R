# MR, 13.07.2018 
########### complementery script: Cervical cancer incidence in Switzerland and type ratios


# copied from R_ce folder: hpvprev_natsal3.R
# install.packages("weight")

require(sciplot)
require(weights)
require(plotrix)
require(plyr)
require(Hmisc)
require(PropCIs)
require(ggplot2)

############# script for HPV prevalences from Natsal. Based on the HPV_prev-1.R script (old)


load("C:/Users/mriesen/Dropbox/Natsal/Natsal-3/natsal3urine.RData")
load("C:/Users/mriesen.CAMPUS/Dropbox/Natsal/Natsal-3/natsal3urine.RData")

natsal3 <- natsal3urine 

agegroups

natsal3$alter <- natsal3$dage
natsal3$agecat[natsal3$alter>=15 & natsal3$alter<=19] <- "15-19"
natsal3$agecat[natsal3$alter>=20 & natsal3$alter<=24] <- "20-24"
natsal3$agecat[natsal3$alter>=25 & natsal3$alter<=29] <- "25-29"
natsal3$agecat[natsal3$alter>=30 & natsal3$alter<=39] <- "30-39"
# natsal3$agecat[natsal3$alter>=40 & natsal3$alter<=64] <- "40-64"
# natsal3$agecat[natsal3$alter>=65 & natsal3$alter<=84] <- "65-84"


## The max age for urine sample is 44 y.

table(natsal3$alter, natsal3$hpv_positive)

summary(natsal3$hpv_positive[natsal3$agecat=="40-64"])
summary(natsal3$hpv_positive[natsal3$agecat=="65-84"])
#dUK_a <- natsal3[natsal3$tot1yr<=900 & natsal3$tot1yr>=0,]
dUK_a <- natsal3
dUK_a <- dUK_a[dUK_a$agecat=="15-19" | dUK_a$agecat=="20-24"
               |dUK_a$agecat=="25-29"| dUK_a$agecat=="30-39",]


#exclude men
dUK_a <- dUK_a[!dUK_a$rsex=="Male" ,  ]


#exclude the sample which were not tested
dUK_hpv <- dUK_a[!dUK_a$hpv_positive=="Not tested / inadequate sample",]
dUK_hpv <- dUK_a[!dUK_a$HPV_HR=="Not tested" ,  ]
dUK_hpv <- dUK_a[!dUK_a$HPV_LR=="Not tested" ,  ]

#exclude those who don't have a urine weight (89)
dUK_hpv <- na.omit(dUK_hpv)


#define positive as 1 and negative as zero
dUK_hpv$HR[dUK_hpv$HPV_HR=="Positive"] <- 1
dUK_hpv$HR[dUK_hpv$HPV_HR=="Negative"] <- 0

dUK_hpv$LR[dUK_hpv$HPV_LR=="Positive"] <- 1
dUK_hpv$LR[dUK_hpv$HPV_LR=="Negative"] <- 0

dUK_hpv$hpvt[dUK_hpv$hpv_positive=="HPV positive"] <- 1
dUK_hpv$hpvt[dUK_hpv$hpv_positive=="HPV negative"] <- 0

dUK_hpv$hpv16[dUK_hpv$HPV_16=="Positive"] <- 1
dUK_hpv$hpv16[dUK_hpv$HPV_16=="Negative"] <- 0

dUK_hpv$hpv18[dUK_hpv$HPV_18=="Positive"] <- 1
dUK_hpv$hpv18[dUK_hpv$HPV_18=="Negative"] <- 0

dUK_hpv$hpv11[dUK_hpv$HPV_11=="Positive"] <- 1
dUK_hpv$hpv11[dUK_hpv$HPV_11=="Negative"] <- 0

dUK_hpv$hpv6[dUK_hpv$HPV_6=="Positive"] <- 1
dUK_hpv$hpv6[dUK_hpv$HPV_6=="Negative"] <- 0


dUK_hpv$chlam<- ifelse(dUK_hpv$ct_posconfirmed=="positive", c(1), c(0))

#urine weights
mean(dUK_hpv$urine_wt)

w <- dUK_hpv$urine_wt/
  mean(dUK_hpv$urine_wt)
mean(w)

dUK_hpv$weights <- w

###################################### For different age group and hpv tests
agegroups <- c( "15-19", "20-24","25-29", "30-39", "40-64", "65-84") 

nag <- 1:3 #age groups 1 to 3
hpvtested <- c("HR", "LR", "hpvt", "hpv16", "hpv18", "hpv11", "hpv6")

nathpvprv <- array(0, dim=c(length(nag), length(hpvtested),3 ), dimnames = list(agegroups[1:length(nag)], hpvtested,
                                                                                c("PointEst", "Lower", "Upper") ) )

names(dUK_hpv)

for( i in 1:length(nag) ){
  for(j in 1:length(hpvtested) ){
    v1 <- as.numeric(t(dUK_hpv[dUK_hpv$agecat==agegroups[ nag[i] ], ][51+j] ) )
    w1 <- (dUK_hpv$urine_wt[dUK_hpv$agecat==agegroups[ nag[i] ] ]/mean(dUK_hpv$urine_wt[dUK_hpv$agecat==agegroups[ nag[i] ] ]) )
    nathpvprv[i,j,] <- binconf( weighted.mean( v1 ,w1 , na.rm=TRUE) * length(v1) ,
                                length(v1 ) )
    # nathpvprv[i,j,] <- ( weighted.mean( v1 ,w1 , na.rm=TRUE)  )
  }}

save(nathpvprv, file="nathpvprv.RData")

round(nathpvprv*100,2)

nathpvprv[,4,1] / nathpvprv[,3,1] 
nathpvprv[,4,1] / nathpvprv[,1,1] 
nathpvprv[,5,1] / nathpvprv[,1,1] 

nathpvprv[,1,1] / nathpvprv[,3,1] 


