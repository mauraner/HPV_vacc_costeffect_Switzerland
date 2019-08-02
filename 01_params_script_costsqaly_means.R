#MR 01.02.19
#Applies mean values of paramsets_CQok to cost and Qalys parameters

disc <-  mean(paramsets_CQok[,which(ptv_cq=="disc")])
Ct[1,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="Ct_1")])
Ct[2,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="Ct_2")])
Ct[3,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="Ct_3")])
Ct[4,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="Ct_4")])
Ct[5,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="Ct_5")])
Cpt    <-  mean(paramsets_CQok[,which(ptv_cq=="Cpt")])  #costs per negative cercival screen
Cfp    <-  mean(paramsets_CQok[,which(ptv_cq=="Cfp")]) #costs of follow-up of false positive
Cme    <-  mean(paramsets_CQok[,which(ptv_cq=="Cme")]) #cost per mild adverse effect (ae)
Cse    <-  mean(paramsets_CQok[,which(ptv_cq=="Cse")]) #cost per severe adverse effect (ae)
#vaccination costs and price diff 
C4v <-   mean(paramsets_CQok[,which(ptv_cq=="C4v")])
C9v <- mean(paramsets_CQok[,which(ptv_cq=="C9v")])
# utilities weights
q[1,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="q_1")])
q[2,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="q_2")])
q[3,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="q_3")])
q[4,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="q_4")])
q[5,1]  <-  mean(paramsets_CQok[,which(ptv_cq=="q_5")])