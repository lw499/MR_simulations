library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(5)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  options(warn=-1)
  library(geepack);library(MASS);library(ResourceSelection);library(ltmle); library(SuperLearner)
  library(dplyr)
  library(data.table)
  setDTthreads(1)
  
  logit <- function(term) {
    return( ifelse(!is.na(term),log(term/(1-term)),NA) )
  }
  
  EXPIT <- function(term) {
    return( ifelse(!is.na(term),exp(term)/(1+exp(term)),NA) )
  }
  
  source("datagen.R")
  set.seed(1129)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
  n <- 1000
  K <- 4
  alpha0=1; alpha1=1; alpha2=-2; alpha3=1;
  beta0=-2; beta1=-2; beta2=-1;
  theta0=-1; theta1=2; theta2=-2; theta3=2; theta4=2;
  kappa0 = 0; kappa1 = 1; kappa2 = 1; kappa3 = -1; 
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(ind, K=K,
            alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3,
            beta0=beta0, beta1=beta1, beta2=beta2,
            theta0=theta0, theta1=theta1, theta2=theta2, theta3=theta3, theta4=theta4,
            kappa0=kappa0, kappa1=kappa1, kappa2=kappa2, kappa3=kappa3)
  })
  
  dffull <- rbindlist(df)
  
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull$lag_A = ifelse(dffull$t0==0,0,dffull$lag_A)
  
  afit = glm(A ~ L1 + L2 + L1*L2, family = binomial(), data = dffull[dffull$lag_A==0,]) #This is from the data generation mechanism

  dffull$pred_obs = predict(afit, newdata = dffull, type="response")
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
      dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)

  dffullwide = dcast(dffull, id ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","Q"))
  
  tmpdata = dffullwide
  #subset data
  tmpdata$Y_1 = ifelse(tmpdata$Y_0==0,0,tmpdata$Y_1)
  tmpdata$Y_2 = ifelse(!is.na(tmpdata$Y_1) & tmpdata$Y_1==0,0,tmpdata$Y_2)
  tmpdata$Y_3 = ifelse(!is.na(tmpdata$Y_2) & tmpdata$Y_2==0,0,tmpdata$Y_3)
  
  tmpdata$Y_4 = tmpdata$Y_3
  tmpdata$Y_3 = tmpdata$Y_2
  tmpdata$Y_2 = tmpdata$Y_1
  tmpdata$Y_1 = tmpdata$Y_0
  
  tmpdata$C_3 = tmpdata$C_2 = tmpdata$C_1 = tmpdata$C_0 = tmpdata$Y_0 = NULL
  
  tmpdata$id = seq(1,n,by=1)
  tmpdata$pi3 <- tmpdata$pi2 <- tmpdata$pi1 <- tmpdata$pi0 <- NA

  tmpdata$pi0 = tmpdata$pred_obs_0
  tmpdata$pi1 = tmpdata$pi0*tmpdata$pred_obs_1
  tmpdata$pi2 = tmpdata$pi1*tmpdata$pred_obs_2
  tmpdata$pi3 = tmpdata$pi2*tmpdata$pred_obs_3

  ##calculate risk
  mean = NULL
  ind = NA;
  #time 1
  tmp = tmpdata[tmpdata$A_0==tmpdata$Q_0,]
  component1 = tmp$Y_1/tmp$pi0;
  comp1 = sum(component1); 
  if(nrow(tmp)>0){param1 = comp1/n} else{param1=NA}
  #time 2
  tmp  = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 & tmpdata$Y_1==1,]
  component2 = tmp$Y_2/tmp$pi1;
  comp2 = sum(component2)
  if(nrow(tmp)>0){param2 = (comp2)/n} else{param2=NA}
  #time 3
  tmp  = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 & tmpdata$A_2==tmpdata$Q_2 & tmpdata$Y_2==1,]
  component3 = tmp$Y_3/tmp$pi2;
  comp3 = sum(component3)
  if(nrow(tmp)>0){param3 = (comp3)/n} else{param3=NA}
  #time 4
  tmp  = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 & tmpdata$A_2==tmpdata$Q_2 & tmpdata$A_3==tmpdata$Q_3 & tmpdata$Y_3==1,]
  component4 = tmp$Y_4/tmp$pi3; 
  #tmp$A_3=ifelse(tmp$Y_3==0, 0, tmp$A_3); tmp$A_2=ifelse(tmp$Y_2==0, 0, tmp$A_2); tmp$A_1=ifelse(tmp$Y_1==0, 0, tmp$A_1); 
  #tmp$pi3 = ifelse(is.na(tmp$pi3),1,tmp$pi3)
  #fit4 = glm(Y_4 ~ 1, family = binomial(), data = tmp, weights = 1/(pi3)); param4 = plogis(summary(fit4)$coef[1,1])
  comp4 = sum(component4)
  if(nrow(tmp)>0){param4 = (comp4)/n} else{param4=NA}

  #param4
  
  ind = ifelse(comp1==0 | comp2==0 | comp3==0 | comp4==0, 1, ind)
  
  myparam = c(param1, param2, param3, param4, ind)
  #mean = rbind(mean,meantmp)
  #myparam = mean
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"IPW_HT.csv")

stopCluster(cl)
