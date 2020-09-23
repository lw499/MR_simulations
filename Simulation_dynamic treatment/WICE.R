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
  #library(reshape2)  #do not use for data frame only
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
  
  
  ##################
  ######time 3######
  ##################
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 &  
                  tmpdata$A_2==tmpdata$Q_2 &  tmpdata$A_3==tmpdata$Q_3,]; 
  y4fit = glm(Y_4 ~ L1_3 + L2_3 + L1_3*L2_3 + Q_3 + L2_3*Q_3, weight=1/pi3, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 &  
                  tmpdata$A_2==tmpdata$Q_2,]; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response"); 
  y4dat$y4pred = ifelse(y4dat$Y_3==0,0,y4dat$y4pred); 
  
  y4fit = glm(y4pred ~ L1_2 + L2_2 + L1_2*L2_2 + Q_2 + L2_2*Q_2, weight=1/pi2, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1,]; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response"); 
  y4dat$y4pred = ifelse(y4dat$Y_2==0,0,y4dat$y4pred); 
  
  y4fit = glm(y4pred ~ L1_1 + L2_1 + L1_1*L2_1 + Q_1 + L2_1*Q_1, weight=1/pi1, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0,]; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response"); 
  y4dat$y4pred = ifelse(y4dat$Y_1==0,0,y4dat$y4pred); 
  
  y4fit = glm(y4pred ~ L1_0 + L2_0 + L1_0*L2_0, weight=1/pi0, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response"); 
  
  meany4tmp = c(mean(y4dat$y4pred))
  
  meany4 = (meany4tmp)
  #meany4
  
  myparam = cbind(meany4)
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"WICE.csv")

stopCluster(cl)