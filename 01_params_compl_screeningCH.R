# MR, 13.07.2018 
########### complementary script: Screening behaviour in Switzerland, based on SHS 2012 dataset


# copied from R_ce folder: CE_screening_behav_SHS_nv.R

require(sciplot)
require(ggplot2)
require(plotrix)
require(plyr)
require(reshape2)

#load SHS 2012 data set
load("O:/Test_R/R/SHS/shs.RData")


#make age groups
agegroups <- c( "15-19", "20-24","25-29", "30-39", "40-64", "65-84")

shs$agegroup1 <-  ifelse(shs$alter>=15 & shs$alter<=19,agegroups[1], 
                         ifelse(shs$alter>=20 & shs$alter<=24, agegroups[2], 
                                ifelse(shs$alter>=25 & shs$alter<=29,agegroups[3],
                                       ifelse(shs$alter>=30 & shs$alter<=39, agegroups[4],
                                              ifelse(shs$alter>=40 & shs$alter<=64, agegroups[5],
                                                     ifelse(shs$alter>=65 & shs$alter<=84, agegroups[6], NA ))))))


table(shs$agegroup1)
table(shs$alter>85)
######################### screening behaviour

#ever had pap screen
table(shs$tkreb03)
#when last pap screen (year)
table(shs$tkreb12c)


table(shs$tkreb03)
table(shs$tkreb03, shs$agegroup1)
table(shs$tkreb12c)
table(shs$tkreb21) #was it in the last 12 months
table(shs$tkreb03,shs$tkreb12c)
table(shs$tkreb21,shs$tkreb12c)


#interview date
table(shs$intjahr)
table(shs$intmonat)
table(shs$inttag)

#date of last pap smear
table(shs$tkreb12a)
table(shs$tkreb12b) #month
table(shs$tkreb12c) #year


table(shs$tkreb21[shs$tkreb12a=="Weiss nicht" & shs$tkreb21=="Ja" ])
table(shs$tkreb21[shs$tkreb12a=="Weiss nicht" & shs$tkreb21=="Nein" ])
#was it in the last year (for those who did no know the date or year 2011 but not knowing the month)
table(shs$tkreb21)

################################# 
# STEP 1: ever had pap smear yes or no, N
# Step 2: was it in the last year, yes or no, N
#           a) year is 2012 or 2011 and months difference is 0 or negative
#           b) year is <2011 or 2011 but months difference is positive
#           c) said "Weiss nicht" to the question when was the last pap smear (tkreb12a) and answered yes or no


table(shs$alter, shs$tkreb03)



####### STEP 1
#10124 women responded yes or no to the question : ever had pap smear
length( shs$tkreb03[shs$tkreb03=="Ja" | shs$tkreb03=="Nein" ] )
# 82.2% had a pap smear at least once in their life (8318 women)
pr1 <- length( shs$tkreb03[shs$tkreb03=="Ja" ] ) / length( shs$tkreb03[shs$tkreb03=="Ja" | shs$tkreb03=="Nein" ] )


####### STEP 2
#new vector, only positive (answered values for months, other are NA)
shs$tkrebsb2 <- ifelse(shs$tkreb12b<0, NA,shs$tkreb12b )


#new vector, if it was less then a year
shs$papsmly <- ifelse(shs$tkreb12c==2012 ,1,
                      ifelse( (shs$tkreb12c==2011 & (shs$tkrebsb2- shs$intmonat) <= 0 )  , 1,
                              ifelse(shs$tkreb12a=="Weiss nicht" & shs$tkreb21=="Ja", 1,
                                     ifelse(  (shs$tkreb12c>1 & shs$tkreb12c<2011 ), 0, 
                                              ifelse( shs$tkreb12c==2011 & (shs$tkrebsb2-shs$intmonat) > 0 ,0,
                                                      ifelse( shs$tkreb12a=="Weiss nicht" & shs$tkreb21=="Nein" , 0,  NA ) ) ) ) ) )


tempd <- shs[shs$tkreb12c>=2011,]
table(tempd$tkrebsb2 - tempd$intmonat )
table(is.na(tempd$tkrebsb2 - tempd$intmonat ) )
table(is.na(tempd$tkrebsb2 ) )
#462 did indicate 2011 but did not indicate the month or whether it was in the last year..
table(shs$papsmly)
# 7881 could say whether it was in the last year or not (462 indicated 2011 but could not say if it was in the last year or not (<12months))
sum(table(shs$papsmly))
# 43% had a papsmear in the last year 
table(shs$papsmly)[2] / sum(table(shs$papsmly))

table(shs$papsmly, shs$tkreb03)

############# Summary table
shs$papsm <- ifelse(shs$tkreb03=="Ja", 1, ifelse(shs$tkreb03=="Nein",0, shs$tkreb03 ))

cdata <- ddply(shs[shs$tkreb03=="Ja" | shs$tkreb03=="Nein", ], .(agegroup1 ), summarise,
               N= length(papsm) ,
               Nly= length(papsmly[!is.na(papsmly)] ),
               mean_overall= mean(papsm, na.rm=TRUE),
               mean_lasty= as.numeric(prop.test(table(papsmly)[2],Nly)[4]),
               mean_oly = mean_overall * mean_lasty,
               lCI_1= prop.test(table(papsm)[2],N)$conf.int[1],
               uCI_1= prop.test(table(papsm)[2],N)$conf.int[2],
               lCI_2= prop.test(table(papsmly)[2],Nly)$conf.int[1],
               uCI_2= prop.test(table(papsmly)[2],Nly)$conf.int[2],
               lCI_3= prop.test(mean_oly*Nly,Nly)$conf.int[1],
               uCI_3= prop.test(mean_oly*Nly,Nly)$conf.int[2] )


cdata
cdata[-6,]

sum(cdata$Nly)

# write.csv(cdata[-6,], file = "sceeningrates_2012SHS_v2.csv")


cdata

####################plot screening rates over age groups

cdata <- cdata[-6,]

plotdir <- "O:/Cost_effectiveness project/Preparation_drafts/figs/"

pdf(paste(plotdir, "screeningSHS", Sys.Date(), ".pdf", sep=""), width=8, height=6)
par(mar = c(5, 4, 1, 1)) 

plot(NA, xlim=c(1, length(cdata[,1] ) ), ylim=c(0,max(cdata[,4:12]))  , 
     ylab= "% screened",
     xlab="",
     frame.plot=FALSE, xaxt="n" , main="screening proportions")

axis(side=1, at=1:length(cdata[,1] ), labels=FALSE)
text(seq(1, length(cdata[,1] ), by=1), par("usr")[3] - 0.005, labels = cdata[,1], srt = 0, pos = 1, xpd = TRUE)

points(1:length(cdata[,1]),cdata[,4], col="black", lwd=2, lty=3, pch=19)
lines(1:length(cdata[,1]),cdata[,4], col="black", lwd=2, lty=3)
segments(x0=1:length(cdata[,1]),
         y0= cdata[,7],
         x1 = 1:length(cdata[,1]),
         y1= cdata[,8], col="black", lwd=2)


points(1:length(cdata[,1])+0.01,cdata[,5], col="black", lwd=2, lty=3, pch=19)
lines(1:length(cdata[,1])+0.01,cdata[,5], col="black", lwd=2, lty=4)
segments(x0=1:length(cdata[,1])+0.01,
         y0= cdata[,9],
         x1 = 1:length(cdata[,1])+0.01,
         y1= cdata[,10], col="black", lwd=2)


points(1:length(cdata[,1]),cdata[,6], col="black", lwd=2, lty=3, pch=17)
lines(1:length(cdata[,1]),cdata[,6], col="black", lwd=2)
segments(x0=1:length(cdata[,1]),
         y0= cdata[,11],
         x1 = 1:length(cdata[,1]),
         y1= cdata[,12], col="black", lwd=2)

legend(1, min(cdata[,3:5], na.rm=TRUE) , legend=c( "ever","last year (ever +)", "last year overall"), 
       pch=c(19,19,17), col="black", lwd=2, lty=c(3,4,1),  bty="n", title="Pap screening" )


mtext('age groups', side = 1, outer = TRUE, line = -2)
# mtext('HPV prevalence', side = 2, outer = TRUE, line = 0)

dev.off()
