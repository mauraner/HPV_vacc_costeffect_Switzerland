# MR, 12.07.2018 
########### Make two sexual activity groups based on Swiss data (01_params_compl_partnerchangerate1_CHsexbehav.R )

# Based on SHS 2012 and Natsal'3 data sets
# For different sex, 2 sexual activity groups and age groups

# copied from R_ce folder: MLE_CE_model_v3.R

library(stats4)
#the shs dataset including the transformation of totnew partners for het. into new partner in last year is included in the
# data set called shs_2, created with the script CE_compl_partnerchangerate_1_retrieveSwissdata.R, stored as : shs_5_inclnewppy.RData

load("tab_params/shs_5_inclnewppy.RData")
summary(shs_2)

# Agegroups were defined before and are :
agegroups <- c( "15-19", "20-24","25-29", "30-39", "40-64", "65-84") #max 72 for natsal3, and 74 for shs
sex <- c("f", "m")

agegroups

summary(shs_2$newppy)
table(shs_2$rsex2)
table(shs_2$agegroup1)
# max(shs_a3_s_f)
# max(shs_a3_s_m)


shs_a <- array(NA, dim=c(length(agegroups), length(sex) ), dimnames = list(agegroups, sex )  )

#in male first age group, there is one indicating 80 partners. I will delete this. I'll also delete one indicatin 99 partners (male, 42-64 agr) 
shs_2 <- shs_2[shs_2$newppy<80,]


for(i in 1:length(sex)){
  for(j in 1:length(agegroups)){
    v1 <- shs_2$newppy[ shs_2$agegroup1==agegroups[j] & shs_2$rsex2==sex[i] ]
    assign(paste0("shs_a", j, "_s_", sex[i],  sep=""),v1 )
    v2 <- (shs_2$wght[ shs_2$agegroup1==agegroups[j] & shs_2$rsex2==sex[i] ]/
             mean(shs_2$wght[ shs_2$agegroup1==agegroups[j] & shs_2$rsex2==sex[i] ] ) )
    assign(paste0("whgt_a", j, "_s_", sex[i] , sep=""), v2 )
  }} #º gives, for each sex and age groups, individual vectors of number of new contacts in the last years and weights

l_hetn1y_f <- list(shs_a1_s_f, shs_a2_s_f, shs_a3_s_f,  shs_a4_s_f,  shs_a5_s_f,  shs_a6_s_f)
l_hetn1y_m <- list(shs_a1_s_m,   shs_a2_s_m,  shs_a3_s_m,   shs_a4_s_m,   shs_a5_s_m,  shs_a6_s_m)

l_w_hetn1y_f <- list( whgt_a1_s_f,whgt_a2_s_f, whgt_a3_s_f,   whgt_a4_s_f,    whgt_a5_s_f,    whgt_a6_s_f )
l_w_hetn1y_m <- list( whgt_a1_s_f,whgt_a2_s_f, whgt_a3_s_m,    whgt_a4_s_m,   whgt_a5_s_m,   whgt_a6_s_m )


big_l <- list(l_hetn1y_f, l_hetn1y_m)
big_w <- list(l_w_hetn1y_f, l_w_hetn1y_m)

big_l[[1]][[3]] #first is sex, second is agegroup

for(j in 1:length(agegroups)){
  for( i in 1:length(sex)){
    print(max(big_l[[i]][[j]]))
  }}

max(big_l[[2]][[3]])    
big_l[[2]][[3]]

#vector with all
big_l_vect <- c()
big_l_vect <- unlist(big_l)
big_w_vect <- unlist(big_w)

mean(big_l[[1]][[2]]    )
##### MLE of two poisson processes #####

############################ First would be to calculates all in one step, but the prop. in low risk group turned out to be too high, 
#here I show only calculations for two steps


#################################### In two steps (first proportions in low activity group, an then the means in low and high for each age groups and sex) 
#first estimate the prop in high or low risk group
v2  <-  c()
f <- function(x,a,m1,m2) {
  a*dpois(x,m1)+(1-a)*dpois(x,m2)
} # function which has three parameter to define two Poisson distributions (a = prop in low risk group, 1-a = proportion in high risk group, m1 = mean of low risk group, m2=mean of high risk group )
nll <- function(a, m1t, m2t) {
  v2 = -sum( big_w_vect *log( f( big_l_vect , a, m1t,  m2t  ) )  ) 
  return( v2  ) }

(est1 <- mle(nll,start=list(a=0.8, m1t=0.5, m2t=5), method="SANN"  )  )

summary(est1)

afix <- est1@coef[1] # the first parameter is set

#second estimate the optimal estimate for each sex and age groups given the proportion calculated above
vm <-  array(NA, dim=c(length(agegroups)-1, length(sex), 2 ), dimnames = list(agegroups[-7], sex, c("m1", "m2") )  )
v3  <-  array(NA, dim=c(length(agegroups)-1, length(sex) ), dimnames = list(agegroups[-7], sex)  )

stv <- 0.5 # is the initial value for m1 

f <- function(x,a,m1,m2) {
  a*dpois(x,m1)+(1-a)*dpois(x,m2)
}
nll <- function( m1a1f, m1a1m, m1a2f, m1a2m, m1a3f, m1a3m, m1a4f,m1a4m, m1a5f, m1a5m, m1a6f, m1a6m, 
                 m2a1f, m2a1m, m2a2f, m2a2m, m2a3f, m2a3m, m2a4f,m2a4m, m2a5f, m2a5m, m2a6f, m2a6m) {
  vm[1,1,1] <- m1a1f ;  vm[1,2,1]  <- m1a1m 
  vm[2,1,1] <- m1a2f ;  vm[2,2,1]  <- m1a2m 
  vm[3,1,1] <- m1a3f ;  vm[3,2,1]  <- m1a3m 
  vm[4,1,1] <- m1a4f ;  vm[4,2,1]  <- m1a4m 
  vm[5,1,1] <- m1a5f ;  vm[5,2,1]  <- m1a5m 
  vm[6,1,1] <- m1a6f ;  vm[6,2,1]  <- m1a6m
  
  vm[1,1,2] <- m2a1f ;  vm[1,2,2]  <- m2a1m 
  vm[2,1,2] <- m2a2f ;  vm[2,2,2]  <- m2a2m 
  vm[3,1,2] <- m2a3f ;  vm[3,2,2]  <- m2a3m 
  vm[4,1,2] <- m2a4f ;  vm[4,2,2]  <- m2a4m 
  vm[5,1,2] <- m2a5f ;  vm[5,2,2]  <- m2a5m 
  vm[6,1,2] <- m2a6f ;  vm[6,2,2]  <- m2a6m 
  
  for(j in 1:(length(agegroups)-1) ){
    for( i in 1:length(sex)){
      v3[j,i] = -sum( big_w[[i]][[j]] *log( f( big_l[[i]][[j]] , afix, vm[j,i,1],   vm[j,i,2]  ) )  ) 
    }}
  return( sum(v3 ) ) }


(est2 <- mle(nll,start=list( m1a1f = stv, m1a1m =stv , m1a2f =stv, m1a2m=stv, m1a3f =stv, m1a3m=stv, m1a4f=stv,m1a4m =stv,
                             m1a5f =stv, m1a5m =stv, m1a6f=stv, m1a6m=stv, m2a1f=5, m2a1m=8,m2a2f=5,m2a2m=8, m2a3f=3, m2a3m =6,
                             m2a4f=3,m2a4m =4, m2a5f=2, m2a5m =2, m2a6f =2, m2a6m=2 ), method="Nelder-Mead" ))


summary(est2)

sbehav <- summary(est2)

est2@coef

sbehav <- as.data.frame( sbehav@coef[,1:2]) 
sprop <- as.data.frame(summary(est1)@coef[,1:2])
sbehav2 <- rbind(sprop[1,],sbehav)

write.csv(sbehav2, file="sbehav_ch3.csv")

### plot
i <- 1
j <- 3
data <- big_l[[i]][[j]]
par(mfrow=c(1,1))
h <- hist(data,breaks=seq(-0.5,max(data)+0.5,1), ylab="frequency", xlab="new heterosexual partners per year", main="Histogram of data")
plot(h$mids,h$density,log="y",xlim=c(0,30),ylim=c(1e-4,1),xlab="new heterosexual partners per year",ylab="frequency",axes=F, 
     main="Maximum Likelihood Estimation (MLE)\n of two Poisson processes")
x <- 0:max(data)
axis(1)
axis(2)
lines(x,f(x, afix, coef(est2)[ ((j*2)-2)+i ], coef(est2)[ ((j*2)-2)+i + 12 ]),lty=1,col="red")
lines(x, afix* dpois(x,coef(est2)[ ((j*2)-2)+i ]),lty=2,col="red")
lines(x,(1-afix)*dpois(x,coef(est2)[((j*2)-2)+i + 12 ]),lty=3,col="red")


## the mle is very dependent on the start values for m2. I will put women and men together to try havine bigger groups and betters estimates
for(j in 1:length(agegroups)){
  v1 <- shs_2$newppy[ shs_2$agegroup1==agegroups[j]  ]
  assign(paste0("shs_b", j,    sep=""),v1 )
  v2 <- (shs_2$wght[ shs_2$agegroup1==agegroups[j] ]/
           mean(shs_2$wght[ shs_2$agegroup1==agegroups[j]] ) )
  assign(paste0("whgt_b", j,  sep=""), v2 )
}

l_hetn1y <- list( shs_b1,  shs_b2,  shs_b3,  shs_b4,  shs_b5,shs_b6 )
l_w_hetn1y <- list( whgt_b1,   shs_b2,    shs_b3,    shs_b4,  shs_b5, shs_b6 )

### PLOT THE PROBABILITY DENSITY
i <-2
data <- l_hetn1y[[i]]
par(mfrow=c(1,1))
h <- hist(data,breaks=seq(-0.5,max(data)+0.5,1), ylab="frequency", xlab="new heterosexual partners per year", main="Histogram of data")
plot(h$mids,h$density,log="y",xlim=c(0,30),ylim=c(1e-4,1),xlab="new heterosexual partners per year",ylab="frequency",axes=F, 
     main="Maximum Likelihood Estimation (MLE)\n of two Poisson processes")
x <- 0:max(data)
axis(1)
axis(2)
lines(x,f(x,coef(est1)[1], coef(est2)[i], coef(est2)[i+6]),lty=1,col="red")
lines(x,coef(est1)[1]*dpois(x,coef(est2)[i]),lty=2,col="red")
lines(x,(1-coef(est1)[1])*dpois(x,coef(est2)[i+6]),lty=3,col="red")

lines(x,f(x,coef(est1)[1], 0.9, 6),lty=1,col="red")
lines(x,coef(est1)[1]*dpois(x,coef(est2)[i]),lty=2,col="red")
lines(x,(1-coef(est1)[1])*dpois(x,coef(est2)[i+6]),lty=3,col="red")




##### for single age group
vm <-  array(NA, dim=c(length(agegroups),2 ), dimnames = list(agegroups,  c("m1", "m2") )  )
v4  <-  c()
j <- 3

f <- function(x,a,m1,m2) {
  a*dpois(x,m1)+(1-a)*dpois(x,m2)
}
nll <- function( m1a1, m2a1) {
  
  v4[j] = -sum( l_w_hetn1y[[j]] *log( f(l_hetn1y[[j]] , afix, m1a1,  m2a1  ) )  ) 
  
  return( sum(v4 ) ) }

(est2 <- mle(nll,start=list( m1a1=1,  m2a1=10), method="Nelder-Mead"  )  )

summary(est2)

sbehav <- summary(est2)





