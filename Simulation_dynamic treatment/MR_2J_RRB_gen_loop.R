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
  library(dplyr);library(nleqslv);
  library(data.table); library(glm2)
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
  
  afit = glm2(A ~ L1 + L2 + L1*L2, family = binomial(), data = dffull[dffull$lag_A==0,]) #This is from the data generation mechanism
  
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
  
  tmpindA = NULL; tmpindY = NULL;
  tmppredA = NULL; tmpindQ = NULL;
  tmpfollow = NULL;
  for(j in 1:(K)) #dwide[,c(paste("overallcen","_",4, sep=""))]
  {
    tmpindA = c(tmpindA, paste("A","_",j-1, sep=""))
    tmpindQ = c(tmpindQ, paste("Q","_",j-1, sep=""))
    tmpfollow = c(tmpfollow, paste("Fol","_",j-1, sep=""))
    tmpindY = c(tmpindY, paste("Y","_",j, sep=""))
    tmppredA = c(tmppredA, paste("pi","_",j-1, sep=""))
    tmpdata[ , tmppredA] <- NA
  }
  
  tmpdata$Fol_0 = as.numeric(tmpdata$A_0==tmpdata$Q_0)
  for(i in 1:(K-1))
  {
    tmppi = select(tmpdata, c(paste("Fol","_",i-1, sep="")))[[1]] * as.numeric(select(tmpdata, c(paste("A","_",i, sep="")))[[1]] == select(tmpdata, c(paste("Q","_",i, sep="")))[[1]])
    set(tmpdata, NULL , paste0("Fol_",i), tmppi)
  }
  ####################
  ####################
  y5dat = tmpdata[ tmpdata[[!is.na(tmpindY[K-1])]] & tmpdata[[tmpindY[K-1]]] == 1 & tmpdata[[tmpfollow[K]]] == 1,] 
  y5dat$Survtmp = select(y5dat, c(paste("Y","_",K, sep="")))[[1]]
  formulatmp = paste0('Survtmp ~', ' L1_', K-1, ' + L2_', K-1, ' + Q_', K-1, ' + L1_', K-1, '*L2_', K-1, ' + L2_', K-1, '*Q_', K-1)
  
  res_list <- vector(mode="list", length=K)
  res_list[[K]] = glm2(formulatmp, family = binomial(), data = y5dat)
  
  y5predlist <- vector(mode="list", length=K)
  y5predlist[[K]] = predict(res_list[[K]], newdata = y5dat)
  
  for(r in (K-1):2)
  {
    y5dat$Survtmp = select(y5dat, c(paste("Y","_",K, sep="")))[[1]]
    pitmp = 1;
    for(j in (K-1):r)
    {
      pitmp = pitmp*select(y5dat, c(paste("pred_obs","_",j, sep="")))[[1]]
    }
    y5dat$pitmp = pitmp
    formulatmp = paste0('Survtmp ~', ' L1_', r-1, ' + L2_', r-1, ' + Q_', r-1, ' + L1_', r-1, '*L2_', r-1, ' + L2_', r-1, '*Q_', r-1)
    res_list[[r]] = glm2(formulatmp, family = binomial(), offset = y5predlist[[r+1]], weight = 1/pitmp, data = y5dat)
    y5predlist[[r]] = predict(res_list[[r]], newdata = y5dat)
  }
  y5dat$Survtmp = select(y5dat, c(paste("Y","_",K, sep="")))[[1]]
  pitmp = 1;
  for(j in (K-1):0)
  {
    pitmp = pitmp*select(y5dat, c(paste("pred_obs","_",j, sep="")))[[1]]
  }
  y5dat$pitmp = pitmp
  formulatmp = paste0('Survtmp ~', ' L1_', 0, ' + L2_', 0, ' + L1_', 0, '*L2_', 0)
  res_list[[1]] = glm2(formulatmp, family = binomial(), offset = y5predlist[[2]], weight = 1/pitmp, data = y5dat)
  y5predlist[[1]] = predict(res_list[[1]], newdata = y5dat)
  
  #### new loops
  for(m in (K-1):3)
  {
  y5predlist_loop <- vector(mode="list", length=m+1);
  y5predlist <- vector(mode="list", length=m+1); y5predlistp <- vector(mode="list", length=m+1)
  y5dat = tmpdata[ tmpdata[[!is.na(tmpindY[m-1])]] & tmpdata[[tmpindY[m-1]]] == 1 & tmpdata[[tmpfollow[m]]] == 1,] 
  for(r in (m+1):1)
  {
    y5predlist[[r]] = predict(res_list[[r]], newdata = y5dat)
    y5predlist_loop[[r]] = y5predlist[[r]]
    y5predlistp[[r]] = predict(res_list[[r]], newdata = y5dat, type="response")
    y5predlistp[[r]] = ifelse(select(y5dat, c(paste("Y","_",m, sep="")))[[1]]==0, 0, y5predlistp[[r]])
  }
  y5dat$Survtmp = y5predlistp[[m]]
  formulatmp = paste0('Survtmp ~', ' L1_', m-1, ' + L2_', m-1, ' + Q_', m-1, ' + L1_', m-1, '*L2_', m-1, ' + L2_', m-1, '*Q_', m-1)
  
  res_listloop <- vector(mode="list", length=m)
  res_listloop[[m]] = glm2(formulatmp, family = binomial(), data = y5dat)
  
  y5predlist_loop <- vector(mode="list", length=m)
  y5predlist_loop[[m]] = predict(res_listloop[[m]], newdata = y5dat)
  
  for(r in (m-1):2)
  {
    y5dat$Survtmp = y5predlistp[[r]]
    pitmp = 1;
    for(j in (m-1):r)
    {
      pitmp = pitmp*select(y5dat, c(paste("pred_obs","_",j, sep="")))[[1]]
    }
    y5dat$pitmp = pitmp
    formulatmp = paste0('Survtmp ~', ' L1_', r-1, ' + L2_', r-1, ' + Q_', r-1, ' + L1_', r-1, '*L2_', r-1, ' + L2_', r-1, '*Q_', r-1)
    res_listloop[[r]] = glm2(formulatmp, family = binomial(), offset = y5predlist_loop[[r+1]], weight = 1/pitmp, data = y5dat)
    y5predlist_loop[[r]] = predict(res_listloop[[r]], newdata = y5dat)
  }
  y5dat$Survtmp = y5predlistp[[1]]
  pitmp = 1;
  for(j in (m-1):0)
  {
    pitmp = pitmp*select(y5dat, c(paste("pred_obs","_",j, sep="")))[[1]]
  }
  y5dat$pitmp = pitmp
  formulatmp = paste0('Survtmp ~', ' L1_', 0, ' + L2_', 0, ' + L1_', 0, '*L2_', 0)
  res_listloop[[1]] = glm2(formulatmp, family = binomial(), offset = y5predlist_loop[[2]], weight = 1/pitmp, data = y5dat)
  y5predlist_loop[[1]] = predict(res_listloop[[1]], newdata = y5dat)
  
  res_list = res_listloop; res_listloop = NULL;
  y5predlist_loop = NULL;
  }
  
  ##second last loop
  m=2
  y5predlist_loop <- vector(mode="list", length=m+1); y5predlistp <- vector(mode="list", length=m+1)
  y5dat = tmpdata[ tmpdata[[!is.na(tmpindY[m-1])]] & tmpdata[[tmpindY[m-1]]] == 1 & tmpdata[[tmpfollow[m]]] == 1,] 
  for(r in (m+1):1)
  {
    y5predlist_loop[[r]] = predict(res_list[[r]], newdata = y5dat)
    y5predlistp[[r]] = predict(res_list[[r]], newdata = y5dat, type="response")
    y5predlistp[[r]] = ifelse(select(y5dat, c(paste("Y","_",m, sep="")))[[1]]==0, 0, y5predlistp[[r]])
  }
  y5dat$Survtmp = y5predlistp[[m]]
  formulatmp = paste0('Survtmp ~', ' L1_', m-1, ' + L2_', m-1, ' + Q_', m-1, ' + L1_', m-1, '*L2_', m-1, ' + L2_', m-1, '*Q_', m-1)
  
  res_listloop <- vector(mode="list", length=m)
  res_listloop[[m]] = glm2(formulatmp, family = binomial(), data = y5dat)
  
  y5predlist_loop <- vector(mode="list", length=m)
  y5predlist_loop[[m]] = predict(res_listloop[[m]], newdata = y5dat)
  
  y5dat$Survtmp = y5predlistp[[1]]
  pitmp = 1;
  for(j in (m-1):0)
  {
    pitmp = pitmp*select(y5dat, c(paste("pred_obs","_",j, sep="")))[[1]]
  }
  y5dat$pitmp = pitmp
  formulatmp = paste0('Survtmp ~', ' L1_', 0, ' + L2_', 0, ' + L1_', 0, '*L2_', 0)
  res_listloop[[1]] = glm2(formulatmp, family = binomial(), offset = y5predlist_loop[[2]], weight = 1/pitmp, data = y5dat)
  y5predlist_loop[[1]] = predict(res_listloop[[1]], newdata = y5dat)
  
  res_list = res_listloop; res_listloop = NULL;
  y5predlist_loop = NULL;
  
  ##last loop
  m=1
  y5predlist_loop <- vector(mode="list", length=m+1); y5predlistp <- vector(mode="list", length=m+1)
  y5dat = tmpdata[tmpdata[[tmpfollow[m]]] == 1,] 
  for(r in (m+1):1)
  {
    y5predlist_loop[[r]] = predict(res_list[[r]], newdata = y5dat)
    y5predlistp[[r]] = predict(res_list[[r]], newdata = y5dat, type="response")
    y5predlistp[[r]] = ifelse(select(y5dat, c(paste("Y","_",m, sep="")))[[1]]==0, 0, y5predlistp[[r]])
  }
  y5dat$Survtmp = y5predlistp[[1]]
  pitmp = 1;
  for(j in (m-1):0)
  {
    pitmp = pitmp*select(y5dat, c(paste("pred_obs","_",j, sep="")))[[1]]
  }
  y5dat$pitmp = pitmp
  res_listloop <- vector(mode="list", length=m)
  #y5predlist_loop <- vector(mode="list", length=m)

  formulatmp = paste0('Survtmp ~', ' L1_', 0, ' + L2_', 0, ' + L1_', 0, '*L2_', 0)
  res_listloop[[1]] = glm2(formulatmp, family = binomial(), weight = 1/pitmp, data = y5dat)
  y5dat = tmpdata
  pred = predict(res_listloop[[1]], newdata = tmpdata, type="response")
  
  meany4 = mean(pred)
  meany4
  
  myparam = cbind(meany4)
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"MR_2J_RRB.csv")

stopCluster(cl)