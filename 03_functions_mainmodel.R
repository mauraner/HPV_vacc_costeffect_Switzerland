# MR, 12.07.2018 
########### Dynamic HPV model based on ordinary differential equations

# Used to model, for each single HPV strain: S-I-R-V + CIN2-CIN3- + CANCER STAGES (undiagnosed and diagnosed) + deaths
# For different sex, sexual activity groups and age groups

# copied from R_ce folder: CE_model_v6.R

####################### Model

#load packages
require(deSolve)
require(rootSolve)

###########################  transmission: sexual mixing between risk groups and age groups (rho_sra):

# Sexual mixing between risk groups AND between age group is here based on assortativity indexes: eps_r and eps_a

rho_sra <- function(eps_r, eps_a) {
  rho <-   array(NA, c(length(sex), length(riskgroups),length(riskgroups), length(agegroups), length(agegroups)), 
                 dimnames = list(sex,riskgroups, paste(riskgroups, "p", sep="_"), agegroups, paste(agegroups, "p", sep="_") ) )
  for (i in 1:length(sex)) {
    for (j in 1:length(riskgroups)){
      for(a in 2:length(agegroups)){
        
        for(jp in 1:length(riskgroups)){
          for(ap in 2:length(agegroups)){
            
            rho[i,j,jp,a,ap] <-   
              (( eps_a*(ifelse(a==ap,1,0))  + (1-eps_a)  
                 *( (sum( C[ifelse(i==1,2,1), 1:length(riskgroups),ap] * N[ifelse(i==1,2,1), 1:length(riskgroups), ap] )) /
                      (sum( C[ifelse(i==1,2,1) ,  1:length(riskgroups), 1:length(agegroups)] 
                            * N[ifelse(i==1,2,1),  1:length(riskgroups), 1:length(agegroups)] ))) )
               *( eps_r*ifelse(j==jp,1,0)  + (1-eps_r)  
                  * ( ( C[ifelse(i==1,2,1), jp, ap] * N[ifelse(i==1,2,1),jp, ap] ) /
                        (sum(C[ifelse(i==1,2,1), 1:length(riskgroups) , ap] 
                             * N[ifelse(i==1,2,1), 1:length(riskgroups) , ap]  ))) ) )
          }}
      }} }
  return(rho)
}


rho <- rho_sra(eps_r, eps_a)

rho[1,1,1,,2]

#adjsted partner change rate:
Cp <- rho

for (i in 1:length(sex)) {
  for (j in 1:length(riskgroups)){
    for(a in 2:length(agegroups)){
      
      for(jp in 1:length(riskgroups)){
        for(ap in 2:length(agegroups)){
          
          Cp[i,j,jp,a,ap] <- C[i,j,a] * sqrt( (C[ifelse(i==1, 2, 1), jp, ap ] * rho[ifelse(i==1, 2, 1), jp,j, ap,a] 
                                               * N[ifelse(i==1, 2, 1), jp,ap] )/ 
                                                (C[i, j, a ] * rho[i, jp,j, a, ap] * N[i, j,a] )  ) 
        }}}}}

Cp[is.na(Cp)] <- 0
rho[is.na(rho)] <- 0

###before running the  model, define:
## the hpv_type choosen:
#                         w =  which(hpvtype=="hpv16")
## disease_outcome choosen:
#                         y =  which(dis_outcome=="cervix") 
lsg <- length(sex)
lrg <- length(riskgroups)
lagr <-  length(agegroups)

HPV_1<- function(t, x, params, verbose=FALSE)
{
  with(as.list(params), {
    dX <- array(0, dim=c(length(comp),length(sex),length(riskgroups),length(agegroups)) )
    ### convert the state vector into an array object
    X <- array(x, dim=c(length(comp),length(sex),length(riskgroups),length(agegroups))) 
    lambda <- array(NA, dim=c(length(sex),length(riskgroups),length(agegroups)) )
    beta <- c(betap[1,w], betap[2,w])
    for (i in 1:lsg){
      ii <- i%%2+1
      for (j in 1:lrg){
        jj <- j%%2+1
        for(a in 1:lagr){
          if(a==1){ 
            aa <- 1
            a2 <- 2
          }else{
            aa <- 0
            a2 <- a}
          
          lambda[i,j,a] <- sum(beta[i]*  rho[i,j,1:lrg,a,1:lagr ]* Cp[i,j,1:lrg,a,1:lagr ]*
                                             ( X[2,ii, 1:lrg, 1:lagr]/
                                                 N[ii, 1:lrg, 1:lagr] ) )
          
          ###1 Susceptible    
          dX[1,i,j,a] <- (  
            ( ( (1-aa) *g[a2-1]*X[1,i,j,a2-1]) + ( aa *(1- pv[i,w,vs]*ve[w] ) * g[length(agegroups)]*N[i,j,length(agegroups) ]  ) )
            - lambda[i,j,a]*X[1,i,j,a] 
            + omega[i,w]*X[12,i,j,a]
            + (1-alpha[w]) * gamma[i,w]* X[2,i,j,a] # percentage that does not seroconvert
            - g[a]*X[1,i,j,a]
            + (s[i,a,y]+ mt*s[i,a,y]) *ssens *X[3,i,j,a] # screened women CIN2
            + (s[i,a,y]+ mt*s[i,a,y]) *ssens* X[4,i,j,a] # screened women CIN3
            + pr2s[w,y]*X[3,i,j,a] #regression from cin2 to normal
            + pr3s[w,y]*X[4,i,j,a] #regression from cin3 to normal
            - m*X[1,i,j,a] + m*rgr[j]*(X[1,i,j,a]+X[1,i,jj,a])
            + sum(pi[1:3,y]*z[1:3,y] *X[c(5,7,9),i,j,a] )
          ) 
          
          ###Infected compartements:
          #2 HPV16
          dX[2,i,j,a] <-  (lambda[i,j,a]* (X[1,i,j,a] + tau[w,y]*X[12,i,j,a]) 
                           - pri1[i, w, y]*X[2,i,j,a]
                           -  pri2[i, w, y]*X[2,i,j,a]
                           - gamma[i,w]*X[2,i,j,a]
                           + (1-aa)* g[a2-1]*X[2,i,j,a2-1]
                           - g[a]*X[2,i,j,a]   
                           - m*X[2,i,j,a] + m*rgr[j]*(X[2,i,j,a]+X[2,i,jj,a]))
          
          #3 CIN2
          dX[3,i,j,a] <-  (pri1[i,w,y]*X[2,i,j,a] 
                           - pr23[w,y]*X[3,i,j,a] 
                           - pr2c[w]*X[3,i,j,a] 
                           - pr2s[w,y]*X[3,i,j,a]
                           + (1-aa)*  g[a2-1]*X[3,i,j,a2-1] # aging in
                           - g[a]*X[3,i,j,a]   #aging out
                           +  pr32*X[4,i,j,a] 
                           - (s[i,a,y] + mt*s[i,a,y])*ssens*X[3,i,j,a] # screened women
                           - m*X[3,i,j,a] + m*rgr[j]*(X[3,i,j,a]+X[3,i,jj,a]))
          
          #4 CIN3
          dX[4,i,j,a] <-  (pri2[i,w,y]*X[2,i,j,a]
                           + pr23[w,y]*X[3,i,j,a] 
                           - pr3c[w]*X[4,i,j,a] 
                           - pr32*X[4,i,j,a]
                           - pr3s[w,y]*X[4,i,j,a]
                           + (1-aa)* g[a2-1]*X[4,i,j,a2-1]
                           - g[a]*X[4,i,j,a] 
                           - (s[i,a,y]+ mt*s[i,a,y])*ssens*X[4,i,j,a] 
                           - m*X[4,i,j,a] + m*rgr[j]*(X[4,i,j,a]+X[4,i,jj,a])) 
          
          #5 undiagnosed local cervical cancer
          dX[5,i,j,a] <-  (pr2c[w]*X[3,i,j,a] #progression CIN2 to ULCC
                           + pr3c[w]*X[4,i,j,a] #progression CIN3 to ULCC
                           + (1-aa)* g[a2-1]*X[5,i,j,a2-1] #aging in
                           - g[a]*X[5,i,j,a]  #aging out
                           - (prlcrc+z[1,y]+um[1]) * X[5,i,j,a] #progression local to regional CC and stage spec. prob of diagn.and mortality
                           - m*X[5,i,j,a] + m*rgr[j]*(X[5,i,j,a]+X[5,i,jj,a])) 
          
          #6 diagnosed local cervical cancer
          dX[6,i,j,a] <-  ( (1-pi[1, y] )* z[1,y]*X[5,i,j,a] #in: stage spec. prob of diagn. (times 1- probability of cure with treatm.)
                            + (1-aa)* g[a2-1]*X[6,i,j,a2-1] #aging in
                            - g[a]*X[6,i,j,a]  #aging out
                            - dm[1] * X[6,i,j,a] #diagnosed mortality
                            - m*X[6,i,j,a] + m*rgr[j]*(X[6,i,j,a]+X[6,i,jj,a])) 
          
          #7 undiagnosed regional cervical cancer
          dX[7,i,j,a] <-  (prlcrc*X[5,i,j,a] 
                           + (1-aa)*g[a2-1]*X[7,i,j,a2-1] #aging in
                           - g[a]*X[7,i,j,a]  #aging out
                           - (prrcdc+z[2,y]+um[2]) * X[7,i,j,a] #progression local to regional CC and stage spec. prob of diagn.and mortality
                           - m*X[7,i,j,a] + m*rgr[j]*(X[7,i,j,a]+X[7,i,jj,a])) 
          
          #8 diagnosed regional cervical cancer
          dX[8,i,j,a] <-  ( (1-pi[2,y] )* z[2,y]*X[7,i,j,a] #in: stage spec. prob of diagn. (times 1- probability of cure with treatm.)
                            + (1-aa)* g[a2-1]*X[8,i,j,a2-1] #aging in
                            - g[a]*X[8,i,j,a]  #aging out
                            - dm[2] * X[8,i,j,a] # diagnosed mortality
                            - m*X[8,i,j,a] + m*rgr[j]*(X[8,i,j,a]+X[8,i,jj,a])) 
          
          #9 undiagnosed distant cervical cancer
          dX[9,i,j,a] <-  (prrcdc*X[7,i,j,a] 
                           + (1-aa)* g[a2-1]*X[9,i,j,a2-1] #aging in
                           - g[a]*X[9,i,j,a]  #aging out
                           - (z[3,y]+um[3]) * X[9,i,j,a] #progression local to regional CC and stage spec. prob of diagn.and mortality
                           - m*X[9,i,j,a] + m*rgr[j]*(X[9,i,j,a]+X[9,i,jj,a])) 
          
          #10 diagnosed distant cervical cancer
          dX[10,i,j,a] <-  ( (1-pi[3,y] )* z[3,y]*X[9,i,j,a] #in: stage spec. prob of diagn. (times 1- probability of cure with treatm.)
                             + (1-aa)* g[a2-1]*X[10,i,j,a2-1] #aging in
                             - g[a]*X[10,i,j,a]  #aging out
                             - dm[3] * X[10,i,j,a] #diagnosed mortality
                             - m*X[10,i,j,a] + m*rgr[j]*(X[10,i,j,a]+X[10,i,jj,a])) 
          
          #11 deaths from CC
          dX[11,i,j,a] <-   ( sum(dm[1:3]* X[c(6,8,10),i,j,a])
                              +sum(um[1:3]* X[c(5,7,9),i,j,a])
                              + (1-aa)* g[a2-1]*X[11,i,j,a2-1]#aging in
                              - g[a]*X[11,i,j,a] )  #aging out)
          
          
          #12 R
          dX[12,i,j,a] <-  (  alpha[w] * gamma[i,w]* X[2,i,j,a]
                              - omega[i,w] * X[12,i,j,a]
                              - lambda[i,j,a]* tau[w,y]*X[12,i,j,a]
                              + (1-aa)* g[a2-1]*X[12,i,j,a2-1]
                              - g[a]*X[12,i,j,a]  
                              - m*X[12,i,j,a] + m*rgr[j]*(X[12,i,j,a]+X[12,i,jj,a]) ) 
          
          #13 Vaccinated 
          dX[13,i,j,a] <- ( ( ( (1-aa)* g[a2-1]*X[13,i,j,a2-1]) + ( (aa)* (pv[i,w,vs]*g[length(agegroups)]*ve[w] *N[i,j,length(agegroups) ] ) ) )
                            - g[a]*X[13,i,j,a] 
                            - m*X[13,i,j,a] + m*rgr[j]*(X[13,i,j,a]+X[13,i,jj,a]) )
          
          # setTxtProgressBar(pba, a)
          # if(verbose){
          # cat("time: ", paste(t, collapse = ";"), '\n')
          # cat("i,j,a: ", paste(i,j,a, collapse = ";"), '\n') }
        }}}
    
    out <- as.vector(dX)
    
    return(list(res = out))
    
  })
}



