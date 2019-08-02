# MR, 13.07.2018 
########### Complementary script. Define age group assortativity, based on Natsal'3 dataset

# copied from R_ce folder: CE_ai_agegroupsassortativtiy_v2.R

require(sciplot)
require(ggplot2)
require(plotrix)
require(plyr)
require(reshape2)

#Load Natsal'3 dataset
# load("C:/Users/mriesen/Dropbox/Maurane/Data/Natsal-3/natsal3.RData")
load("C:/Users/mriesen.CAMPUS/Dropbox/Maurane/Data/Natsal-3/natsal3.RData")


summary(natsal3$rsex)
summary(natsal3$dage)
summary(natsal3$tot1yr)
hist(natsal3$tot1yr)

agegroups <- c("15-19", "20-24","25-29", "30-39", "40-64", "65-84", "85+") #max 72 for natsal3, and 74 for shs
natsal3$agegroup1 <-  ifelse(natsal3$dage>=15 & natsal3$dage<=19,agegroups[1], 
                             ifelse(natsal3$dage>=20 & natsal3$dage<=24, agegroups[2], 
                                    ifelse(natsal3$dage>=25 & natsal3$dage<=29,agegroups[3],
                                           ifelse(natsal3$dage>=30 & natsal3$dage<=39, agegroups[4],
                                                  ifelse(natsal3$dage>=40 & natsal3$dage<=64, agegroups[5],
                                                         ifelse(natsal3$dage>=65 & natsal3$dage<=84, agegroups[6],
                                                                ifelse(natsal3$dage>=85, agegroups[7], NA   ))))))) 

natsal3$agegroup1 <- as.factor(natsal3$agegroup1)
summary(natsal3$agegroup1)

# we will look for assorativity between 18-24, 25-31, 32-41 and 42+
names(natsal3)

#r1patage: age on 1st occasion with most recent partner (from partner)
#rafsmr respondent's age in completed years at 1st sex with most recent partner

table(natsal3$r1ptage)
table(natsal3$rafsmr)

hist( natsal3$rafsmr - natsal3$r1ptage  , breaks=100)

table(natsal3$agegroup1)
summary(natsal3$agegroup1)

# subdata which includes response for r1ptage and rafsmr
natsal3_s1 <- natsal3[natsal3$r1ptage <995 & natsal3$r1ptage >-1 &natsal3$rafsmr <98 & natsal3$rafsmr >-1,  ]

# Group age of partner and age of 1st sex with the partner into age groups
natsal3_s1$agegroup_pmr <- ifelse(natsal3_s1$r1ptage>=15 & natsal3_s1$r1ptage<=19,agegroups[1], 
                                  ifelse(natsal3_s1$r1ptage>=20 & natsal3_s1$r1ptage<=24, agegroups[2], 
                                         ifelse(natsal3_s1$r1ptage>=25 & natsal3_s1$r1ptage<=29,agegroups[3],
                                                ifelse(natsal3_s1$r1ptage>=30 & natsal3_s1$r1ptage<=39, agegroups[4],
                                                       ifelse(natsal3_s1$r1ptage>=40 & natsal3_s1$r1ptage<=64, agegroups[5],
                                                              ifelse(natsal3_s1$r1ptage>=65 & natsal3_s1$r1ptage<=84, agegroups[6],
                                                                     ifelse(natsal3_s1$r1ptage>=85, agegroups[7], NA   ))))))) 


natsal3_s1$agegroup_rmr <-  ifelse(natsal3_s1$rafsmr>=15 & natsal3_s1$rafsmr<=19,agegroups[1], 
                                   ifelse(natsal3_s1$rafsmr>=20 & natsal3_s1$rafsmr<=24, agegroups[2], 
                                          ifelse(natsal3_s1$rafsmr>=25 & natsal3_s1$rafsmr<=29,agegroups[3],
                                                 ifelse(natsal3_s1$rafsmr>=30 & natsal3_s1$rafsmr<=39, agegroups[4],
                                                        ifelse(natsal3_s1$rafsmr>=40 & natsal3_s1$rafsmr<=64, agegroups[5],
                                                               ifelse(natsal3_s1$rafsmr>=65 & natsal3_s1$rafsmr<=84, agegroups[6],
                                                                      ifelse(natsal3_s1$rafsmr>=85, agegroups[7], NA   ))))))) 

table(natsal3_s1$agegroup_rmr)
summary(natsal3_s1$r1ptage)

table(natsal3_s1$r1ptage>15)

table(natsal3_s1$agegroup_pmr)


library("circlize")
require("migest")
require("plyr")
require("reshape2")

############################## Contact made by a group and recieved
dt1 <- as.data.frame( cbind(natsal3_s1$agegroup_rmr,natsal3_s1$agegroup_pmr ) )
names(dt1) <- c("respondent", "partner")

dt2 <- as.data.frame( table(dt1$respondent,dt1$partner) )
chordDiagram(x = dt2)

col1 <- colorRampPalette(c("grey20", "lightgrey"))(6)
col3 <- colorRampPalette(c("darkblue", "lightblue"))(6)

dt3 <- dt2
dt3$cat <- ifelse( dt2$Var1 == dt2$Var2 , "same", "diff" ) 


#define colors of links
g.col.l <- ifelse( dt3$cat=="same" & dt3$Var1==agegroups[1] , col3[1], 
                   ifelse(dt3$cat=="same" & dt3$Var1==agegroups[2] , col3[2],
                          ifelse(dt3$cat=="same" & dt3$Var1==agegroups[3] , col3[3],
                                 ifelse(dt3$cat=="same" & dt3$Var1==agegroups[4] , col3[4],
                                        ifelse(dt3$cat=="same" & dt3$Var1==agegroups[5] ,  col3[5],
                                               ifelse(dt3$cat=="same" & dt3$Var1==agegroups[1] , col3[6],
                                                      "lightgrey" ) )))))

#firt row plot with entering and exiting contacts (not grouped)
#pdf("Migrationcircos_agemix_v4.pdf", width=6, height=6)

chordDiagram(dt2, grid.col= col3, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  

circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)
#dev.off()

#take the assortativity
dt4 <- dcast(dt2, Var1~Var2, value = 'Freq') 
row.names(dt4) <- dt4[,1]
dt4 <- dt4[,2:7]
dt5 <- dt4
i <- 1
for(i in 1:6){
  dt5[,i] <- dt4[,i] / colSums(dt4)[i]
}

dt5 <- as.matrix(dt5)

mean( diag(dt5) )

#this is not a symetrical matrix. I will sum the contacts between agegroups
################ First I take only the contacts made:
dt6 <- dt4 


dt6[lower.tri(dt6, diag=FALSE)] <- 0
dt6 <- as.matrix(dt6)

dt7 <- melt(dt6)

#plot of only contacts made
pdf("Migrationcircos_agemix_v4sym_contactexit.pdf", width=6, height=6)

chordDiagram(dt7, grid.col= col3, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  
circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)

dev.off()

dt8 <- dt6

tempd <- t(dt8)

dt8[lower.tri(dt8)] <- tempd[lower.tri(tempd)]
dt9 <- dt8

isSymmetric.matrix(dt9)
isSymmetric.matrix(dt8)

i <- 1
for(i in 1:6){
  dt9[,i] <- dt8[,i] / colSums(dt8)[i]
}

dt9 <- as.matrix(dt9)
colSums(dt9)
mean(  diag(dt9) )

#weighted for number of contacts per age cat: 
colSums(dt8)

weighted.mean(diag(dt9), colSums(dt8))


################ Second I sum up the contacts made and double the diagonals:
dt6 <- dt4 
# diag(dt4) <- 0

tempm <- (t(dt4))

symtri <- dt4[upper.tri(dt4)] + tempm[upper.tri(tempm) ] 
dt6[upper.tri(dt6)] <- symtri
dt6[lower.tri(dt6)]  <- 0
dt6 <- as.matrix(dt6)
diag(dt6) <- 2*diag(dt6)

dt7 <- melt(dt6)

# plot of contacts made and recieved
pdf("Migrationcircos_agemix_v5sym_exitenter.pdf", width=6, height=6)

chordDiagram(dt7, grid.col= col3, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  


circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)

dev.off()

dt8 <- dt6

tempd <- t(dt8)

dt8[lower.tri(dt8)] <- tempd[lower.tri(tempd)]
dt9 <- dt8

isSymmetric.matrix(dt9)
isSymmetric.matrix(dt8)

i <- 1
for(i in 1:6){
  dt9[,i] <- dt8[,i] / colSums(dt8)[i]
}

dt9 <- as.matrix(dt9)
colSums(dt9)

mean(  diag(dt9) )

#weighted for number of contacts per age cat: 
colSums(dt8)

weighted.mean(diag(dt9), colSums(dt8))
# 

diag(dt8)/rowSums(dt8)


