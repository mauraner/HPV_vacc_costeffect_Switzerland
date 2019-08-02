# MR,  nv 23.05.19
########### complementery script: Cervical cancer incidence in Switzerland and type ratios

# copied from R_ce folder: CE_cancer_incidence.R
dirl <- "tab_params/"

# retrieve incidence data from National Institut for Cancer Epidemiology and Registration (NICER)
# http://www.nicer.org/en/statistics-atlas/cancer-incidence/

#contains data for whole Switzerland, men and women, different cancer, number of cases (observed/estimated),
# Crude rate per 100.000 person-year, rate per 100.000 person-year standardised by age (EU standard) and LL95%CI and UL%CI
require("readr")

CH_incid <- read_csv("tab_params/Swiss_cancer_incid.csv",   locale = locale(encoding = "UTF8"))
CH_incid <- as.data.frame(CH_incid)

summary(CH_incid)

#C53 is the cervical cancer incidence
table(CH_incid$`ICD-10`=="C53")

dat_cc <- CH_incid[CH_incid$`ICD-10`=="C53",]

# plot crude cervical cancer incidence
plot(dat_cc$Période, dat_cc$`Taux brut`, ylim=c(0,12), xlim=c(1990,2014) , xlab="year", ylab="crude incidence rate / 100'000 person-year", 
     main="Crude cervical cancer incidence in Switzerland", frame=F, pch=16, cex=1.5, xaxt="n")
axis(1, at= seq(1990,2014, by=2) )

pdf("cc_crude_CH.pdf", width=8, height=6)
plot(dat_cc$Période, dat_cc$`Taux brut`, ylim=c(0,12), xlim=c(1990,2014) , xlab="year", ylab="incidence rate / 100'000 person-year", 
     main="Crude cervical cancer incidence in Switzerland", frame=F, pch=16, cex=1, xaxt="n")
axis(1, at= seq(1990,2014, by=2) )
lines(dat_cc$Période, dat_cc$`Taux brut`)

lines(x=2004:2014, y=rep(mean( dat_cc$`Taux brut`[dat_cc$Période>2003]), 11), lwd=2, col="darkblue" )
text(2011,mean( dat_cc$`Taux brut`[dat_cc$Période>2003])+3, 
     paste("mean incidence, \nfrom 2004-2014: \n", round(mean( dat_cc$`Taux brut`[dat_cc$Période>2003]),2), sep="" ), col="darkblue" )

dev.off()


# plot standardised cervical cancer incidence, over years (ESP 1967)
# pdf("cc_standardised_CH_vthesis.pdf", width=8, height=6)
plot(dat_cc$Période, dat_cc$`Taux standardisé`, ylim=c(0,12), xlim=c(1990,2014) , xlab="year", ylab="incidence rate / 100'000 person-year", 
     main="age-standardised cervical cancer incidence in Switzerland", frame=F, pch=16, cex=1, xaxt="n")
axis(1, at= seq(1990,2014, by=2) )
lines(dat_cc$Période, dat_cc$`Taux standardisé`)

segments(x0=dat_cc$Période,
         y0= dat_cc$`lower 95%`,
         x1 = dat_cc$Période,
         y1= dat_cc$`upper IC 95%`, col="black")
lines(x=2004:2014, y=rep(mean( dat_cc$`Taux standardisé`[dat_cc$Période>2003]), 11), lwd=2, col="darkblue" )
text(2011,mean( dat_cc$`Taux standardisé`[dat_cc$Période>2003])+3, 
     paste("mean incidence, \nfrom 2004-2014: \n", round(mean( dat_cc$`Taux standardisé`[dat_cc$Période>2003]),2), sep="" ), col="darkblue" )
# dev.off()

#numer of cases 2004-2014 
mean( dat_cc$`Nombre de cas`[dat_cc$Période>2003] )



### It seems to go down naturally until 2004 and then reaches like a plateau from then onwards. To fit my data, I will take a 10 year
# period from 2004-2014
totinc_cc <- mean( dat_cc$`Taux brut`[dat_cc$Période>2003])

# Retrieve a SD which includes the min and max estimates, centered around the mean estimate 
min(dat_cc$`Taux brut`[dat_cc$Période>2003])
max(dat_cc$`Taux brut`[dat_cc$Période>2003])

minB <- (totinc_cc*2) - max(dat_cc$`upper IC 95%`[dat_cc$Période>2003])
maxB <- max(dat_cc$`upper IC 95%`[dat_cc$Période>2003])

CHccSD <- totinc_cc -minB

#Per HPV types (based on Clifford et al. 2003 BJC systematic review, region Europe, including 33 studies and 3336 women)
# I take the sum of only high risk hpy types: 16,18,31,33,35,39,45,51,52,56,58,59,68,73,82. First I take those included in the model,
#then I add the other HR types and 100- the propostion with any hpv type to that
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")


hpvtypprop <- hpvtype
names(hpvtypprop) <- hpvtype

#other HR in the study, grouped together
otherrestinghr <- c(0.1,0.6,0.4,0.0,0.3) #types 59,56,68,82,73,

# hpvtypprop <- c(56.0, 17.5, 4.2, 4.4, 0.5, mean(otherresting), 2.9, 0.0, 0.5, 0.6, mean(otherresting),
#                 0.1, 0.1, mean(otherresting),0, 0.0)

#looking at % of multiple infections:
names(hpvtypprop) <- hpvtype
hpvtypprop <- c(56.0, 17.5, 4.2, 4.4, 0.5, 0.3, 2.9, 0.0, 0.5, 0.8, 0.0,0.0, sum(otherrestinghr) )

sum(hpvtypprop) + (100-85.9) #you'd have 2.6% of co-infection cases in the cervical cancer

hpvtypprop <- hpvtypprop/ ( sum(hpvtypprop) ) # I do not take into account co-infections and assume that 100% of cervical cancer cases are caused by HPV


cc_inc_pertype2 <- hpvtypprop * totinc_cc

sum(hpvtypprop)

names(cc_inc_pertype2) <- hpvtype

sum(cc_inc_pertype2)

# write.csv(cc_inc_pertype2, file = "cc_inc_pertype4e6.csv")


############### Now with a more recent source (2010 de Sanjosé, Lance Oncol.)
hpvtype <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv39", "hpv45",
             "hpv51", "hpv52", "hpv58", "hpv6", "hpv11", "other_HR")
hpvtypprop <- hpvtype
names(hpvtypprop) <- hpvtype

notmodelledhpvtypesHR <- c(3, 5, 1,1, 44,32,15,2,3,13,4,1,1,10)
names(notmodelledhpvtypesHR) <- c("26","30", "34", "39","44","56", "59", "66","67", "68", "68b", "70", "91" , "underter")

names(hpvtypprop) <- hpvtype
hpvtyppropnb <- c(1348, 150, 69, 117, 46, 27, 80, 28, 40 , 27,  0.0,0.0, sum(notmodelledhpvtypesHR) )

hpvtypprop <- hpvtyppropnb / sum(hpvtyppropnb)

cc_inc_pertype1 <- hpvtypprop * totinc_cc

sum( cc_inc_pertype1[c(1,2,3,4,7,9,10,11,12)] ) / sum(cc_inc_pertype1)
sum(hpvtyppropnb)

write.csv(cc_inc_pertype1, file = "tab_params/cc_inc_pertype4e6.csv")

sum( cc_inc_pertype1[c(5,6,8,13)] ) # is the sum of the non vaccine type remaining

sum( cc_inc_pertype1[-c(1,2)] ) 

#total 2058 patients, excl: hpv6, 11, 42, 51, 53, 66, 73. total 2067 infection types

############################# Further look at cervical cancer incidence in Switzerland: looking at differences between language regions
CH_incid_inclang <- read_csv("tab_params/Swiss_cancer_incid_incl_langreg.csv",
                             locale = locale(encoding = "UTF-8"))
CH_incid_inclang <- as.data.frame(CH_incid_inclang)

names(CH_incid_inclang)

dat_cc_inclang <- CH_incid_inclang[CH_incid_inclang$`ICD-10`=="C53",]

# pdf("cc_inc_langugageregion_withregline.pdf", width=8, height=6)
plot(dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse alémanique"], dat_cc_inclang$`Taux standardisé`[dat_cc_inclang[,1]=="Suisse alémanique"], 
     ylim=c(0,12), xlim=c(1990,2014) , xlab="year", ylab="incidence rate / 100'000 person-year", 
     main="Standardised cervical cancer incidence in Switzerland", frame=F, pch=16, cex=1.5, xaxt="n", col="grey")
axis(1, at= seq(1990,2014, by=2) )

points(dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse franco-italienne"]+0.2, dat_cc_inclang$`Taux standardisé`[dat_cc_inclang[,1]=="Suisse franco-italienne"], 
       pch=16, cex=1.5, xaxt="n", col="darkblue")

legend(1990,4, legend=c("German", "French and Italian"), pch=16, col=c("grey", "darkblue"), bty="n", title="language region:" )

segments(x0=dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse alémanique"],
         y0= dat_cc_inclang$`lower 95%`[dat_cc_inclang[,1]=="Suisse alémanique"],
         x1 = dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse alémanique"],
         y1= dat_cc_inclang$`upper IC 95%`[dat_cc_inclang[,1]=="Suisse alémanique"], col="grey")
segments(x0=dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse franco-italienne"]+0.2,
         y0= dat_cc_inclang$`lower 95%`[dat_cc_inclang[,1]=="Suisse franco-italienne"],
         x1 = dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse franco-italienne"]+0.2,
         y1= dat_cc_inclang$`upper IC 95%`[dat_cc_inclang[,1]=="Suisse franco-italienne"], col="darkblue")

abline(m1, col="grey")
abline(m2, col="darkblue")


# dev.off()

m1 <- glm(dat_cc_inclang$`Taux standardisé`[dat_cc_inclang[,1]=="Suisse alémanique"] ~dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse alémanique"])
m2 <- glm(dat_cc_inclang$`Taux standardisé`[dat_cc_inclang[,1]=="Suisse franco-italienne"] ~dat_cc_inclang$Période[dat_cc_inclang[,1]=="Suisse franco-italienne"])


#################### Other cancer stages incidences linked with cervical cancer, type ratio from different sources

################################ CIN 3 type ratios
#Guan et al 2012, hpv type ratio for all hpv positive samples (not specific data for Europe)
typeratiocin3 <-       c(6259, 791, 337, 941, 869, 1215, 987, 341, 299, 203, 509, 234, 149)
typeratiocin3tested <- c(10753, 10669, 9348, 10296, 9635, 10386, 9652, 9481, 9105, 8977, 9195, 9390, 7406)
GUAN_typeratiocin3 <- typeratiocin3 /typeratiocin3tested
names(GUAN_typeratiocin3) <- paste( "hpv", c(16,18,45,33,58,31,52,35,39,59,51,56,68), sep="")
sum(GUAN_typeratiocin3 )
#reorder hpv types according to model order and group other HR together
nGUAN_typeratiocin3 <-  GUAN_typeratiocin3[match( hpvtype, names(GUAN_typeratiocin3) )]
names(nGUAN_typeratiocin3) <- hpvtype
nGUAN_typeratiocin3[11:13] <- c(0,0,sum(GUAN_typeratiocin3[names(GUAN_typeratiocin3) %in% hpvtype ==FALSE] ))
sum(nGUAN_typeratiocin3)
write.csv(nGUAN_typeratiocin3, file = "nGUAN_typeratiocin3.csv")

# Swiss data from the CIN3+plus study, personal communication from Dianne Egli
swisscin3plus <- c(398, 36, 94, 47, 18, 7, 12, 25, 43, 7, 33, 6 )
Swissratiocin3 <- swisscin3plus /sum(swisscin3plus)
names(Swissratiocin3) <- paste( "hpv", c(16,18,31,33,35,39,45,51,52,56,58,59), sep="")
nSwissratiocin3 <-  Swissratiocin3[match( hpvtype, names(Swissratiocin3) )]
names(nSwissratiocin3) <- hpvtype
nSwissratiocin3[11:13] <- c(0,0,mean(Swissratiocin3[names(Swissratiocin3) %in% hpvtype ==FALSE] ))
write.csv(nSwissratiocin3, file = "nSwissratiocin3.csv")

################################ Normal cytology type ratios
#I take estimates for world, 
# typerationormal <- c(6767, 2768, 1462, 1476, 1959, 2490, 2386, 985, 1521,932,1994,1535,869)
# typerationormalpost <- c(33154,32964, 30359, 31121, 30961, 31130, 29660, 29279, 29260, 29163, 28259, 29485, 28192)
# typerationormal / typerationormalpost
# # GUAN_typerationormal <- typerationormal / sum(typerationormal)
# GUAN_typerationormal <- typerationormal / typerationormalpost
# names(GUAN_typerationormal) <- paste( "hpv", c(16,18,45,33,58,31,52,35,39,59,51,56,68), sep="")
# nGUAN_typerationormal <-  GUAN_typerationormal[match( hpvtype, names(GUAN_typerationormal) )]
# names(nGUAN_typerationormal) <- hpvtype
# nGUAN_typerationormal[11:13] <- c(0,0,mean(GUAN_typerationormal[names(GUAN_typerationormal) %in% hpvtype ==FALSE] ) )
# write.csv(nGUAN_typerationormal, file = "nGUAN_typerationormal.csv")


#only Europe, for women with normal cytology
typerationormal <-        c(3344,  1286,  805,    726,   735,  1586, 1019, 476, #europe exept hpv 39,51 and OHR (59, 56, 68) hereafter
                            1521, 1994, 932, 1535, 869) 
typerationormaltested <-  c(14636, 14636, 13434, 13842, 13436, 13877, 12474, 12264, 
                            29260, 28259, 29163, 29485, 28192)

hpvtn <- paste( "hpv", c(  16,      18,    45,     33,    58,    31,   52,  35, 
                           39,    51,  59,   56,   68), sep="")

GUAN_typerationormal <- typerationormal / typerationormaltested
sum(GUAN_typerationormal )

names(GUAN_typerationormal) <- paste( "hpv", c(16,18,45,33,58,31,52,35, 39, 51, 59, 56, 68), sep="")
nGUAN_typerationormal <-  GUAN_typerationormal[match( hpvtype, names(GUAN_typerationormal) )]
names(nGUAN_typerationormal) <- hpvtype

######### other version instead of assuming that up to 1 is other high risk is to account for the prevalence in 59,56,68 reported in GUan 2012 as OHR
nGUAN_typerationormal[11:13] <- c(0,0,sum(GUAN_typerationormal[names(GUAN_typerationormal) %in% hpvtype ==FALSE] ) )
# sum(nGUAN_typerationormal)
nGUAN_typerationormal <- nGUAN_typerationormal/sum(nGUAN_typerationormal)
write.csv(nGUAN_typerationormal, file = "nGUAN_typerationormal.csv")



## version 1, used before, assuming that the rest up to 100% are other HR
# ohr <- 1-sum(GUAN_typerationormal )
# nGUAN_typerationormal[11:13] <- c(0,0,sum(GUAN_typerationormal[names(GUAN_typerationormal) %in% hpvtype ==FALSE] )+ohr )
# write.csv(nGUAN_typerationormal, file = "nGUAN_typerationormal.csv")



