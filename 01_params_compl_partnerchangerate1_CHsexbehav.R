# MR, 12.07.2018 
########### Partner change rates for Switzerland

# Based on SHS 2012 and Natsal'3 data sets
# For different sex, 2 sexual activity groups and age groups


# Transform Swiss SHS data about total number of partner in partner change rate, based on natsal'3 information'

require(sciplot)
require(ggplot2)
require(plotrix)
require(plyr)
require(reshape2)

####################### NATSAL'3
#load dataset
# load("C:/Users/mriesen/Dropbox/Maurane/Natsal/Natsal-3/natsal3.RData")
load("C:/Users/mriesen.CAMPUS/Dropbox/Maurane/Data/Natsal-3/natsal3.RData")

#summary of values of iterest
summary(natsal3$rsex)
summary(natsal3$dage)
summary(natsal3$tot1yr)
hist(natsal3$tot1yr)
summary(natsal3$het1yr)
hist(natsal3$het1yr)

#define age groups wanted
agegroups <- c("15-19", "20-24","25-29", "30-39", "40-64", "65-84", "85+") #max 72 for natsal3, and 74 for shs

natsal3$agegroup1 <- ifelse(natsal3$dage>=15 & natsal3$dage<=19,agegroups[1], 
                            ifelse(natsal3$dage>=20 & natsal3$dage<=24, agegroups[2], 
                                   ifelse(natsal3$dage>=25 & natsal3$dage<=29,agegroups[3],
                                          ifelse(natsal3$dage>=30 & natsal3$dage<=39, agegroups[4],
                                                 ifelse(natsal3$dage>=40 & natsal3$dage<=64, agegroups[5],
                                                        ifelse(natsal3$dage>=65 & natsal3$dage<=84, agegroups[6],
                                                               ifelse(natsal3$dage>=85, agegroups[7], NA   ))))))) 




natsal3$agegroup1 <- as.factor(natsal3$agegroup1)
summary(natsal3$agegroup1)
############################################################# EXPLORING DATA #################
# For age and gender, the total number of partner in the last year, UK
ctot_UK <- ddply(natsal3[!natsal3$het1yr>800,], .(agegroup1, rsex ), summarise,
                 N= length(!is.na(het1yr)) , 
                 mean= mean (het1yr,  na.rm=TRUE), 
                 w.mean= weighted.mean(het1yr,  total_wt, na.rm=TRUE), 
                 sd= sd(het1yr,  na.rm=TRUE),
                 se= sd/ sqrt(N))

# For age and gender, the number of new partners in the last year, UK
cnew_UK <- ddply(natsal3[!natsal3$hetnonew>800,], .(agegroup1, rsex ), summarise,
                 N= length(!is.na(hetnonew)) , 
                 mean= mean (hetnonew,  na.rm=TRUE), 
                 w.mean= weighted.mean(hetnonew,  total_wt, na.rm=TRUE), 
                 sd= sd(hetnonew,  na.rm=TRUE),
                 se= sd/ sqrt(N))


ctot_UK
cnew_UK

newd <- ctot_UK[,1:2]
newd$prop <-  ctot_UK[,5] -cnew_UK[,5] 

newd
newd1 <- cbind(ctot_UK[,c(1,2,5,7) ] , cnew_UK[,c(3,5,7) ] ) #I take the N of tot new (which has usually 1-2 persons less)
newd1 <- newd1[,c(1,2,5,3,4,6,7)]
colnames(newd1) <- c("agegroup1", "sex", "N", "tot_w.mean", "tot_se", "new_w.mean", "new_se")

# New variables with the difference between total and new nubmer of partners
natsal3$diffnewtot <- natsal3$het1yr - natsal3$hetnonew

# number of person for which the differences betwen total and new is 0,1,2, etc. until max possible
dift <- table(natsal3$diffnewtot[natsal3$diffnewtot >= 0 & natsal3$diffnewtot < 990] )
dift / sum(dift)
# number of person for which the differences betwen total and new is 0,1,2 and the rest is grouped together
dift2 <- c(dift[1:3], sum(dift[4:28]) )
dift2/ sum(dift2)

# number of person for which the differences betwen total and new is 0,1,2, etc. until max possible for each value of differences
dift3 <- table(natsal3$diffnewtot[natsal3$diffnewtot >= 0 & natsal3$diffnewtot < 990 &
                                    natsal3$het1yr >= 0 & natsal3$het1yr < 990],
               natsal3$het1yr[natsal3$diffnewtot >= 0 & natsal3$diffnewtot < 990 &
                                natsal3$het1yr >= 0 & natsal3$het1yr < 990]  )

dift4 <- apply( dift3, 1, function(x) c(x[1:6], sum(x[7:length(x)]) ) )
rownames(dift4) <- c(0:5, ">=6")
dift5 <- apply( dift4, 1, function(x) c(x[1:6], sum(x[7:length(x)]) ) )
rownames(dift5) <- c(0:5, ">=6")
#in dift5, columns are number of total partners and rows are the difference between total and new partners
dift6 <- apply( dift5, 2, function(x) x/sum(x)  )

############################################################# FINIHED EXPLORING DATA, START #################

## I need to find out from how many persons I have to make -1, -2, -3 etc. (in relationships) to get the Nnew partners
sex <- c("f", "m")
natsal3$rsex2  <- NA
natsal3$rsex2[natsal3$rsex==1] = "m"
natsal3$rsex2[natsal3$rsex==2] = "f"
table(natsal3$rsex2)
table(natsal3$rsex)

natsal3$diffnewtot <- natsal3$het1yr - natsal3$hetnonew


dtf1 <- natsal3[natsal3$diffnewtot >= 0 & natsal3$diffnewtot < 990 &
                  natsal3$het1yr >= 0 & natsal3$het1yr < 990 ,]


# p_age1 <- array(NA, dim=c(length(agegroups), length(sex),  max(sort(unique(dtf1$het1yr) ) )+1 ,length(unique(dtf1$het1yr)) ) , 
#                 dimnames = list(agegroups, sex, 0:max(sort(unique(dtf1$het1yr) ) ), sort(unique(dtf1$het1yr) ) )  ) 

p_age1 <- array(NA, dim=c(length(agegroups), length(sex),  max(sort(unique(dtf1$het1yr) ) )+1 ,max(sort(unique(dtf1$het1yr) ) )+1  ) , 
                dimnames = list(agegroups, sex, 0:max(sort(unique(dtf1$het1yr) ) ),  0:max(sort(unique(dtf1$het1yr) ) ) )  ) 

for(i in 1:length(sex) ){
  for(j in 1:length(agegroups)){
    dt3 <- table(dtf1$diffnewtot[ dtf1$agegroup1==agegroups[j] & dtf1$rsex2==sex[i]],
                 dtf1$het1yr[ dtf1$agegroup1==agegroups[j] & dtf1$rsex2==sex[i] ])
    dift6 <- apply( dt3, 2, function(x) x/sum(x)  )
    p_age1[j,i, which(sort(unique(dtf1$diffnewtot) ) %in% as.numeric(row.names(dift6)) ),
           which(sort(unique(dtf1$het1yr) ) %in% as.numeric(colnames(dift6)) )] <- dift6
    
  }}

p_age1[3,1,,]
p_age1[4,1,,]
p_age1[j,i,,]

p_age1[1,,,]

##
##################################### Now apply this to the swiss SHS data
require(sciplot)
require(ggplot2)
require(plotrix)
require(plyr)
require(reshape2)

load("O:/Test_R/R/SHS/shs.RData")

summary(shs$sex)
summary(shs$alter)

# SHS: question TAIDS01 asks if the person ever had sex and TAIDS15 asks if the participant had sex in the last 12 months, if yes, 
# it asks how often and in TAIDS17 it asks with how many persons. I will combine taids 01 15 and 17
table(shs$taids01)
table(shs$taids15)
table(shs$taids17)

# shs$nbsexc <- ifelse(shs$taids01=="Nein", "virgin", ifelse(shs$taids15=="Nein", "0", shs$taids17 ) )
shs$nbsexc <- ifelse(shs$taids01=="Nein", 0, ifelse(shs$taids15=="Nein", 0, shs$taids17 ) )

#now exclude homosexual partnerships (if women: if during life they had answered having sex mostly with women but at least,
# with one men, and/or only with women)
shs$orient <- ifelse(shs$sex=="Frau" & shs$taids03=="Haupts. mit Frauen, aber mindestens mit 1 Mann", "hom",
                     ifelse(shs$sex=="Frau" & shs$taids03=="Nur mit Frauen", "hom",
                            ifelse( shs$sex=="Mann" & shs$taids04=="Haupts. mit Männern, aber mindestens mit 1 Frau","hom",
                                    ifelse( shs$sex=="Mann" & shs$taids04=="Nur mit Männern" , "hom" , "het"))))

#make age groups
shs$agegroup1 <-  ifelse(shs$alter>=15 & shs$alter<=19,agegroups[1], 
                         ifelse(shs$alter>=20 & shs$alter<=24, agegroups[2], 
                                ifelse(shs$alter>=25 & shs$alter<=29,agegroups[3],
                                       ifelse(shs$alter>=30 & shs$alter<=39, agegroups[4],
                                              ifelse(shs$alter>=40 & shs$alter<=64, agegroups[5],
                                                     ifelse(shs$alter>=65 & shs$alter<=84, agegroups[6],
                                                            ifelse(shs$alter>=85 , agegroups[7], NA )))))))
table(shs$agegroup1)
table(shs$nbsexc, shs$alter)
summary(shs$nbsexc)
summary(shs$wght)

table(shs$orient)
table(shs$nbsexc)
table(shs$nbsexc<0)
table(shs$nbsexc, shs$alter)

length(shs$nbsexc[!shs$nbsexc<0 & shs$orient=="het"])
# New subset
shs_1 <- shs[!shs$nbsexc<0,]
shs_1 <- shs_1[! shs_1$orien=="hom",]

shs_1$rsex2  <- NA
shs_1$rsex2[shs_1$sex=="Mann"] = "m"
shs_1$rsex2[shs_1$sex=="Frau"] = "f"
table(shs_1$rsex2)
table(shs_1$sex)

table(shs_1$nbsexc, shs_1$agegroup1, shs_1$rsex2)


shs_1 <- shs_1[sample(nrow(shs_1)),]
shs_1$newppy <- NA

p_age1[j,i,,]
p_age2 <- p_age1

#if there is no information about the specific number of tot partner in Natsal3, it refers to the lower max indicated value
for(i in 1:length(sex) ){
  for(j in 1:length(agegroups)){
    tempd <- shs_1$nbsexc[shs_1$rsex2==sex[i] & shs_1$agegroup1==agegroups[j]] 
    for(k in  sort(unique(tempd))[-1] ) {
      closnna <- (which(!is.na(p_age1[j,i,1,]))[which.min( abs( ((as.vector(which(!is.na(p_age1[j,i,1,]) ) )) -1)
                                                                - k ) ) ] )
      pr1 <-   p_age1[j,i,1:closnna, closnna  ]
      #probabtiliy vector replacing NA by zeros
      pr1[is.na(pr1)] <- 0
      shs_1$newppy[shs_1$rsex2==sex[i] 
                   & shs_1$agegroup1==agegroups[j] 
                   & shs_1$nbsexc== k ] <- sample((closnna-1):0, length(tempd[tempd==k]), prob= pr1, replace=T) 
    }         
  }}

# it seems that nbsexc are transformed to NA (instead of staying 0)
summary(shs_1$newppy)
shs_1$newppy[is.na(shs_1$newppy)] <- 0

summary(shs_1$newppy)
summary(shs_1$nbsexc)

length(shs_1$nbsexc)
length(shs_1$newppy)

shs_1$nbsexc[which(is.na(shs_1$newppy) ) ]

shs_2 <- shs_1[,c( which(names(shs_1)=="agegroup1" | names(shs_1)=="newppy" | names(shs_1)=="wght" | names(shs_1)=="rsex2" ) )]

# save(shs_2, file="shs_5_inclnewppy.RData")


##########################
cdata <- ddply(shs_1, .(agegroup1, sex ), summarise,
               N= length(nbsexc) , 
               #tot.mean= mean (nbsexc,  na.rm=TRUE), 
               tot.w.mean = weighted.mean(nbsexc, wght, na.rm=TRUE),
               tot_sd= sd(nbsexc,  na.rm=TRUE),
               tot_se= tot_sd/ sqrt(N),
               #new.mean= mean (newppy,  na.rm=TRUE), 
               new.w.mean = weighted.mean(newppy, wght, na.rm=TRUE),
               new_sd= sd(newppy,  na.rm=TRUE),
               new_se= new_sd/ sqrt(N)
)
cdata


# From the "exploring data" in the beginning
colnames(newd1) <- c("agegroup1", "sex","N" ,"tot_w.mean", "tot_se", "new_w.mean", "new_se")
newd1$data <- "Natsal'3"
newd1$sex[newd1$sex==2] <- "f"
newd1$sex[newd1$sex==1] <- "m"

cdata <- cdata[,c(-5 ,-8)]
colnames(cdata) <- c("agegroup1", "sex","N", "tot_w.mean", "tot_se", "new_w.mean", "new_se")
cdata$sex <- as.character(cdata$sex)
cdata$data <- "SHS based"
cdata$sex[cdata$sex=="Frau"] <- "f"
cdata$sex[cdata$sex=="Mann"] <- "m"

comb_data <- rbind(newd1,cdata )

############################################### PLOT the two data sets (UK and CH)
#colors:
cswissmw <- c("red", "darkred")
cukmw <- c("lightblue", "darkblue")

cswissmw <- c("darkseagreen4", "steelblue4"  )
cukmw <- c("darkseagreen3", "steelblue3"  )

plotdir <- "O:/Cost_effectiveness project/Preparation_drafts/figs/"

pdf( paste(plotdir, "totandnewpartnerperyear_CHUK.pdf", sep= "" ), width=12, height=6)

par(mfrow=c(1,2))
# par(mar=c(7,5,3,1))
plot(NA,  ylim=c(0, max(comb_data$tot_w.mean + 1.96*comb_data$tot_se)+0.8 ), xlim=c(1,length(levels(comb_data$agegroup1) )+1.5  ), 
     xlab= "age groups",
     ylab="mean number of partners",
     cex.lab=1.5, frame.plot=FALSE, xaxt = "n")
axis(1, at=1:6, labels=levels(comb_data$agegroup1))

x_ax1 <- c(1,1.1,2,2.1,3,3.1,4,4.1,5,5.1, 6, 6.1)
x_ax2 <- c(1,2,3,4,5,6)
# Calculted new numbers of partner for Swiss population (from SHS and Natsal'3)
segments(x0=x_ax1,
         y0= comb_data$new_w.mean[comb_data$data=="SHS based"] - 1.96*comb_data$new_se[comb_data$data=="SHS based"],
         x1 = x_ax1,
         y1=  comb_data$new_w.mean[comb_data$data=="SHS based"] + 1.96*comb_data$new_se[comb_data$data=="SHS based"], col="black", lty=1)

points(x_ax1,comb_data$new_w.mean[comb_data$data=="SHS based"], col=cswissmw, pch=19)
for(i in 1:2){
  lines(if(i==1) {x_ax2} else x_ax2+0.1,comb_data$new_w.mean[comb_data$data=="SHS based"& comb_data$sex==c("m", "f")[i]], col=cswissmw[i],  lwd=2)
}
# Total number of partner, SHS 2012
segments(x0=x_ax1,
         y0= comb_data$tot_w.mean[comb_data$data=="SHS based"] - 1.96*comb_data$tot_se[comb_data$data=="SHS based"],
         x1 = x_ax1,
         y1=  comb_data$tot_w.mean[comb_data$data=="SHS based"] + 1.96*comb_data$tot_se[comb_data$data=="SHS based"], col="black", lty=1)
for(i in 1:2){
  lines(if(i==1) {x_ax2} else x_ax2+0.1,comb_data$tot_w.mean[comb_data$data=="SHS based"& comb_data$sex==c("m", "f")[i]], col=cswissmw[i], lwd=2, lty=3)
}
points(x_ax1,comb_data$tot_w.mean[comb_data$data=="SHS based"], col=cswissmw, pch=21, bg="white" )

legend(5,3, legend=c("men","women"), pch=1, lty=3, col=cswissmw,
       bty='n', cex=0.8, title="total last year (SHS 2012)")

legend(5,2.2, legend=c("men","women"), pch=19, lty=1, col=cswissmw,
       bty='n', cex=0.8, title="new last year (calculated)")

# par(mar=c(7,5,3,1))
plot(NA,  ylim=c(0, max(comb_data$tot_w.mean + 1.96*comb_data$tot_se)+0.8 ), xlim=c(1,length(levels(comb_data$agegroup1) )+1.5  ), 
     xlab= "age groups",
     ylab="mean number of partners",
     cex.lab=1.5, frame.plot=FALSE, xaxt = "n")
axis(1, at=1:6, labels=levels(comb_data$agegroup1))

x_ax1 <- c(1,1.1,2,2.1,3,3.1,4,4.1,5,5.1, 6, 6.1)
x_ax2 <- c(1,2,3,4,5,6)

# New numbers of partner from Natsal'3
segments(x0=x_ax1,
         y0= comb_data$new_w.mean[comb_data$data=="Natsal'3"] - 1.96*comb_data$new_se[comb_data$data=="Natsal'3"],
         x1 = x_ax1,
         y1=  comb_data$new_w.mean[comb_data$data=="Natsal'3"] + 1.96*comb_data$new_se[comb_data$data=="Natsal'3"], col="black", lty=1)

points(x_ax1,comb_data$new_w.mean[comb_data$data=="Natsal'3"], col=cukmw, pch=19 )
for(i in 1:2){
  lines(if(i==1) {x_ax2} else x_ax2+0.1,comb_data$new_w.mean[comb_data$data=="Natsal'3"& comb_data$sex==c("m", "f")[i]], col=cukmw[i],  lwd=2)
}
# Total number of partner, Natsal'3
segments(x0=x_ax1,
         y0= comb_data$tot_w.mean[comb_data$data=="Natsal'3"] - 1.96*comb_data$tot_se[comb_data$data=="Natsal'3"],
         x1 = x_ax1,
         y1=  comb_data$tot_w.mean[comb_data$data=="Natsal'3"] + 1.96*comb_data$tot_se[comb_data$data=="Natsal'3"], col="black", lty=1)
for(i in 1:2){
  lines(if(i==1) {x_ax2} else x_ax2+0.1,comb_data$tot_w.mean[comb_data$data=="Natsal'3"& comb_data$sex==c("m", "f")[i]], col=cukmw[i],  lwd=2, lty=3)
}
points(x_ax1,comb_data$tot_w.mean[comb_data$data=="Natsal'3"], col=cukmw, pch=21, bg="white"  )

legend(5,3, legend=c("men","women"), pch=1, lty=3, col=cukmw,
       bty='n', cex=0.8, title="total last year (Natsal-3)")

legend(5,2.2, legend=c("men","women"), pch=19, lty=1, col=cukmw,
       bty='n', cex=0.8, title="new last year (Natsal-3)")


dev.off()

# from the SIRS data set, we know that for youg girls aged 18-24, mean number of new partners per year is : > mean(data_ch)
# [1] 0.5224793, so what we have here is maybe underestimated, the sample size is also smaller



