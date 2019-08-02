## Sensitivity parameter sets to be added on the natural history parameters that are working
# MR, 24.01.19

lsens <-  nsens_ok #nb of combinations for each strat

ptv_cq <- c("disc", "Ct_1", "Ct_2",  "Ct_3", "Ct_4", "Ct_5", "Cpt", "Cfp" , "Cme", "Cse", "C4v", "C9v", 
            "q_1" ,"q_2", "q_3" ,"q_4", "q_5")

paramsets_CQok <- array(0, dim=c(lsens, length(ptv_cq)),  dimnames = list( c(1:lsens),ptv_cq ) ) 


source("01_params_script_costsqaly_baselineparam.R")
#discounting rates
# paramsets_CQok[,which(ptv_cq=="disc")] <- runif(lsens, min=0.0, max=0.06) 
paramsets_CQok[,which(ptv_cq=="disc")] <- 0.03 #does not vary for the multivariate sensitivity analysis
#disease treatement costs
discf <- 0.5
paramsets_CQok[,which(ptv_cq=="Ct_1")] <- rlnorm( lsens,log(Ct[1,1]^2 / sqrt((Ct[1,1]*discf)^2 + Ct[1,1]^2)) , sqrt(log(1 + ((Ct[1,1]*discf)^2 / Ct[1,1]^2))) )
paramsets_CQok[,which(ptv_cq=="Ct_2")] <- rlnorm( lsens,log(Ct[2,1]^2 / sqrt((Ct[2,1]*discf)^2 + Ct[2,1]^2)) , sqrt(log(1 + ((Ct[2,1]*discf)^2 / Ct[2,1]^2))) ) 
paramsets_CQok[,which(ptv_cq=="Ct_3")] <- rlnorm( lsens,log(Ct[3,1]^2 / sqrt((Ct[3,1]*discf)^2 + Ct[3,1]^2)) , sqrt(log(1 + ((Ct[3,1]*discf)^2 / Ct[3,1]^2))) ) 
paramsets_CQok[,which(ptv_cq=="Ct_4")] <- rlnorm( lsens,log(Ct[4,1]^2 / sqrt((Ct[4,1]*discf)^2 + Ct[4,1]^2)) , sqrt(log(1 + ((Ct[4,1]*discf)^2 / Ct[4,1]^2))) ) 
paramsets_CQok[,which(ptv_cq=="Ct_5")] <- rlnorm( lsens,log(Ct[5,1]^2 / sqrt((Ct[5,1]*discf)^2 + Ct[5,1]^2)) , sqrt(log(1 + ((Ct[5,1]*discf)^2 / Ct[5,1]^2))) ) 

paramsets_CQok[,which(ptv_cq=="Cpt")] <- rlnorm( lsens,log(Cpt^2 / sqrt((Cpt*discf)^2 + Cpt^2)) ,
                                                 sqrt(log(1 + ((Cpt*discf)^2 / Cpt^2))) )  #costs per negative cercival screen
paramsets_CQok[,which(ptv_cq=="Cfp")] <- rlnorm( lsens,log(Cfp^2 / sqrt((Cfp*discf)^2 + Cfp^2)) ,
                                                 sqrt(log(1 + ((Cfp*discf)^2 / Cfp^2))) )  #costs of follow-up of false positive
paramsets_CQok[,which(ptv_cq=="Cme")] <- rlnorm( lsens,log(Cme^2 / sqrt((Cme*discf)^2 + Cme^2)) ,
                                                 sqrt(log(1 + ((Cme*discf)^2 / Cme^2))) )  #cost per mild adverse effect (ae)
paramsets_CQok[,which(ptv_cq=="Cse")] <- rlnorm( lsens,log(Cse^2 / sqrt((Cse*discf)^2 + Cse^2)) ,
                                                 sqrt(log(1 + ((Cse*discf)^2 / Cse^2))) )  #cost per severe adverse effect (ae)
#vaccination costs and price diff 
vacpf <- 0.2 #factor for sd of vaccine price
vacpdf <- 10 #factor for sd of vaccine price difference
pricediff <- C9v - C4v
# pricediff <- rlnorm( 1,log(pricediff^2 / sqrt((pricediff*vacpdf)^2 + pricediff^2)) , sqrt(log(1 + ((pricediff*vacpdf)^2 / pricediff^2))) )
# pricediff <- rgamma( 1, shape=pricediff, scale=2 )
pricediff <- runif( lsens,  pricediff*0.1, pricediff*10 )

paramsets_CQok[,which(ptv_cq=="C4v")] <-  rlnorm( lsens,log(C4v^2 / sqrt((C4v*vacpf)^2 + C4v^2)) , sqrt(log(1 + ((C4v*vacpf)^2 / C4v^2))) )
paramsets_CQok[,which(ptv_cq=="C9v")] <- paramsets_CQok[,which(ptv_cq=="C4v")] + pricediff

# C9v <-  rlnorm( 1,log(C9v^2 / sqrt((C9v*vacpf)^2 + C9v^2)) , sqrt(log(1 + ((C9v*vacpf)^2 / C9v^2))) )
# utilities weights
utwf1 <- 0.05
utwf2 <- 0.1

paramsets_CQok[,which(ptv_cq=="q_1")] <- runif(lsens, q[1,1]- q[1,1]*utwf1, q[1,1]+ q[1,1]*utwf1) 
paramsets_CQok[,which(ptv_cq=="q_2")] <- runif(lsens, q[2,1]- q[2,1]*utwf1, q[2,1]+ q[2,1]*utwf1) 

paramsets_CQok[,which(ptv_cq=="q_3")] <- runif(lsens, q[3,1]- q[3,1]*utwf2, q[3,1]+ q[3,1]*utwf2) 
paramsets_CQok[,which(ptv_cq=="q_4")] <- runif(lsens, q[4,1]- q[4,1]*utwf2, q[4,1]+ q[4,1]*utwf2)  
paramsets_CQok[,which(ptv_cq=="q_5")] <- runif(lsens, q[5,1]- q[5,1]*utwf2, q[5,1]+ q[5,1]*utwf2) 


dirl <- "tab_params/" #where to load parameter data from
date_save <- gsub("-", "_", Sys.Date()) 

save(paramsets_CQok, file=paste(dirl,  "paramsets_CQok" , date_save, ".RData", sep="") )


summary(paramsets_CQok[,which(ptv_cq=="C4v")] )
plot(hist(paramsets_CQok[,which(ptv_cq=="q_3")]))
summary(paramsets_CQok[,which(ptv_cq=="q_3")])
summary(paramsets_CQok[,which(ptv_cq=="Ct_1")])
plot(hist(paramsets_CQok[,which(ptv_cq=="Ct_1")]))
summary(paramsets_CQok[,which(ptv_cq=="disc")])

summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="C4v")] )

summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="C9v")] -paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="C4v")] )

summary(paramsets_CQok[,which(ptv_cq=="q_5")] )
summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="q_1")] )
summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="q_2")] )
summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="q_3")] )
summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="q_4")] )
summary(paramsets_CQok[unique(mCEA_sens2b_sub50$Var3),which(ptv_cq=="q_5")] )

summary(rbeta(lsens, estBetaParams(q[1,1],q[1,1]*utwf1)[1], estBetaParams(q[1,1],q[1,1]*utwf1)[2],ncp = 0) )

dim(paramsets_CQok)
summary(paramsets_CQok[,which(ptv_cq=="disc")])

mCEA_sens2b_sub50
