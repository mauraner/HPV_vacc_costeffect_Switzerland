# 01_params_script_workingparamsets.R
# MR, May 2019

# Script to continue the sensitivity analysis with the parameter sets that achieve the calibration
# run the sensitivity analysis on the different scenarios

#directories
dirl <- "tab_params/" #where to load parameter data from
loaddir <- "sens_step1/190528/p" #where the UBE data was stored

source("01_params_script_baselineparam.R")

load( paste(loaddir,7, ".RData", sep=""))

source(paste(loaddir,7, ".RData", sep=""))

load( paste(loaddir,66, ".RData", sep=""))

tosave[[1]]
tosave[[2]]

#####
load(paste(dirl, "paramsets2019_01_26.RData", sep="") ) 
dim(paramsets)

nruns <- 2000 #number of parameter sets

betaused <- array(NA, dim=c(length(hpvtype), nruns), dimnames = list(hpvtype, 1:nruns) )  
pr3cused <- array(NA, dim=c(length(hpvtype), nruns), dimnames = list(hpvtype, 1:nruns) )  

for(i in 1:nruns){
  if( file.exists(paste(loaddir, i,".RData", sep="") ) ) {
    load( paste(loaddir, i,".RData", sep="") )
    betaused[,i] <- unlist(tosave[1])
    pr3cused[,i] <- unlist(tosave[2])
    next
  }}  


summary(betaused[1,])
round( summary(betaused[1,]),2)
summary(pr3cused[1,])

for(i in 1:length(hpvtype)){
  print(hpvtype[i])
  print(summary(betaused[i,]))
  print(summary(pr3cused[i,]))
}

which(is.na(pr3cused[5,]))!= which(is.na(pr3cused[1,]))

table(is.na(pr3cused[5,]))
table(is.na(pr3cused[1,]))

ics <- which(is.na(pr3cused[5,])!= (is.na(pr3cused[1,]))) #parameter sets which did not work for pr3c for hpv 35

# write.table(ics, "ics.txt", sep="\t") 


length(which(is.na(betaused[5,])))
length(which(is.na(pr3cused[5,])))

betaused[5,1087]
pr3cused[5,1087]

#make table with calibrated parameters for each HPV type
tab_calibparams <- as.data.frame( matrix(ncol=2, nrow=length(hpvtype) ) ) 
row.names(tab_calibparams) <- hpvtype
row.names(tab_calibparams)[13] <- "other"
colnames(tab_calibparams) <- c( "$\\beta$ (mean, IQR)", "$\\phi_{3,c}$ (mean, IQR)")
for(i in 1:length(hpvtype) ){
  tab_calibparams[i,1] <- paste( round(summary(betaused[i,]),3)[4], " (", 
                                 round(summary(betaused[i,]),3)[2], "-",
                                 round(summary(betaused[i,]),3)[5], ")", sep="")
  
  tab_calibparams[i,2] <- paste( round(summary(pr3cused[i,]),3)[4], " (", 
                                 round(summary(pr3cused[i,]),3)[2], "-",
                                 round(summary(pr3cused[i,]),3)[5], ")", sep="")
}
tab_calibparams <- tab_calibparams[-c(11,12),]
dirlts <- "C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/"
write.table(tab_calibparams, paste(dirlts,"P_calib_", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)

isnap_1 <- list()
isnap_2 <- list()


for(i in 1:length(hpvtype)){
  isnap_1[[i]] <- (which(is.na(pr3cused[i,])  ) )
  isnap_2[[1]] <- (which(is.na(pr3cused[1,])  ) )
  print(i)
  print(any( (isnap_1[[i]]==isnap_2[[1]]) ))
}

#which nruns are NAs:
for(i in 1:length(hpvtype)){
  print(length(isnap_1[[i]]))
  print(isnap_1[[i]]%in%isnap_1[[1]])
  
}

isnap_1[[5]]
table(isnap_1[[5]]%in%isnap_1[[1]])

for(i in 2:length(hpvtype)){
  print(table( isnap_1[[i]]%in%isnap_1[[1]]) )
}

###############################
unlist(isnap_1[1]) # all not working combination for hpv 16
unlist(isnap_1[5])[1] # one additionnal not working combination

nwps <- c(unlist(isnap_1[5]) )
# nwps <- c(unlist(isnap_1[1]))


# isnap_3 <- (which(is.na(pr3cused[5,])  ) ) #NA including the 34 non working pr3c
# nwps <- c(isnap_3)
length(nwps)
############# New parameter set without the not working combinations

# load(paste(dirl, "paramsets2019_01_15.RData", sep="") ) 
dim(paramsets)

# beta_ok <- array(NA, dim=c(length(hpvtype), nruns-length(nwps) ), dimnames = list(hpvtype, 1: (nruns-length(nwps) ) ) )  
paramsets_ok <- paramsets[1:1116,,]
beta_ok  <- betaused[,1:1116] 
pr3c_ok  <- pr3cused[,1:1116] 

paramsets_ok <- paramsets[-nwps,,]
beta_ok <- betaused[,-nwps] 
pr3c_ok <- pr3cused[,-nwps] 

workingparamset <- (1:2000)[-nwps]
length(workingparamset)

date_save <- gsub("-", "_", Sys.Date()) 

save(paramsets_ok, file=paste(dirl,  "paramsets_ok" , date_save, ".RData", sep="") )
save(beta_ok, file=paste(dirl,  "beta_ok" , date_save, ".RData", sep="") )
save(pr3c_ok, file=paste(dirl,  "pr3c_ok" , date_save, ".RData", sep="") )
save(workingparamset, file=paste(dirl,  "workingparamset" , date_save, ".RData", sep="") )

#table with prior and posterior parameters for each HPV type
# show 19 parameter varied in the 1st step of the sensitivity analysis
sensparam_rownames <- c("$\\alpha$", "$\\gamma_1$", "$\\gamma_2$", "$\\omega$", "$\\nu$", "$\\phi_{12}$", "$\\phi_{13}$",
                        "$\\phi_{23}$", "$\\varphi_{2}$", "$\\varphi_{3}$", "$\\phi_{45}$", "$\\phi_{56}$", "$\\vartheta_{4}$" , "$\\vartheta_{2}$",
                        "$\\vartheta_{6}$", "$\\psi_{4}$", "$\\psi_{6}$", "$\\psi_{8}$", "$\\mu_{4}$", "$\\mu_{6}$", "$\\mu_{8}$",
                        "$\\mu_{5}$","$\\mu_{7}$","$\\mu_{9}$", "$\\sigma_{1}$", "$\\sigma_{2}$", "$\\sigma_{3}$", "$\\sigma_{4}$",
                        "$\\sigma_{5}$", "$\\sigma_{6}$", "$\\sigma_{7}$", "$\\varsigma$", "$\\tau$", "$\\eta$")

length(sensparam_rownames)
# length(ptv)

lhpt <- list(c(1:3),c(4:6), c(7:9), c(10,13))

for(iii in 1:4){
  hpt <- unlist(lhpt[iii])
  tab_sensparams <- as.data.frame( matrix(ncol=length(hpt), nrow= length(ptv) ) ) 
  row.names(tab_sensparams) <- sensparam_rownames
  colnames(tab_sensparams) <- hpvtype[hpt]
  
  for(ii in 1:length(hpt) ){
    for(j in 1:length(ptv)){
      i <- hpt[ii]
      tab_sensparams[j,ii] <- paste( round(summary(paramsets[,j,i]),3)[4], " (", 
                                     round(summary(paramsets[,j,i]),3)[2], "-",
                                     round(summary(paramsets[,j,i]),3)[5], ")", 
                                     "/ ",round(summary(paramsets_ok[,j,i]),3)[4], " (", 
                                     round(summary(paramsets_ok[,j,i]),3)[2], "-",
                                     round(summary(paramsets_ok[,j,i]),3)[5], ")"  ,sep="")
    }}
  
  tab_sensparams <- tab_sensparams[-c(25,26),]
  
  dirlts <- "C:/Users/mriesen.CAMPUS/Documents/R_enroute/Cost_effectiveness project/Preparation_drafts/tab/"
  write.table(tab_sensparams, paste(dirlts,"P_paramsens",iii,"_", Sys.Date(), ".txt" ,sep=""), quote=FALSE, eol="\\\\\n", sep=" & ", col.names = F)
}

