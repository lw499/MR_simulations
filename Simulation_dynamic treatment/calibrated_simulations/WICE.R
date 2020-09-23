library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(10)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  options(warn=-1)
  library(geepack);library(MASS);library(ResourceSelection);library(ltmle); library(SuperLearner)
  library(dplyr)
  library(data.table)
  library(nleqslv)
  #library(reshape2)  #do not use for data frame only
  setDTthreads(1)
  
  logit <- function(term) {
    return( ifelse(!is.na(term),log(term/(1-term)),NA) )
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
  
  
  #testdata = dffull[dffull$lag_A==0,]
  #testdata$logittrans = logit(testdata$pred_obs)
  #testdata$blah = predict(afit, newdata = testdata)
  #testdata$blah = ifelse(testdata$A==1,testdata$blah,-testdata$blah)
  
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
  y4fit1 = glm(Y_4 ~ L1_3 + L2_3 + L1_3*L2_3 + Q_3 + L2_3*Q_3, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 &  
                    tmpdata$A_2==tmpdata$Q_2,]; 
  y4dat$y4pred = predict(y4fit1, newdata = y4dat, type="response"); 
  y4dat$y4pred = ifelse(y4dat$Y_3==0,0,y4dat$y4pred); 
  
  y4fit2 = glm(y4pred ~ L1_2 + L2_2 + L1_2*L2_2 + Q_2 + L2_2*Q_2, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1,]; 
  y4dat$y4pred = predict(y4fit2, newdata = y4dat, type="response"); 
  y4dat$y4pred = ifelse(y4dat$Y_2==0,0,y4dat$y4pred); 
  
  y4fit3 = glm(y4pred ~ L1_1 + L2_1 + L1_1*L2_1 + Q_1 + L2_1*Q_1, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0,]; 
  y4dat$y4pred = predict(y4fit3, newdata = y4dat, type="response"); 
  y4dat$y4pred = ifelse(y4dat$Y_1==0,0,y4dat$y4pred); 
  
  y4fit4 = glm(y4pred ~ L1_0 + L2_0 + L1_0*L2_0, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata; 
  y4dat$y4pred0 = predict(y4fit4, newdata = y4dat, type="response"); 
  y4dat$y4pred1 = predict(y4fit3, newdata = y4dat, type="response"); 
  y4dat$y4pred2 = predict(y4fit2, newdata = y4dat, type="response"); 
  y4dat$y4pred3 = predict(y4fit1, newdata = y4dat, type="response"); 
  
  alldata = y4dat;
  alldata$h3 = logit(alldata$pred_obs_3)
  alldata$h2 = logit(alldata$pred_obs_2)
  alldata$h1 = logit(alldata$pred_obs_1)
  alldata$h0 = logit(alldata$pred_obs_0)
  
  mytmpdat = alldata
  model=function(beta){ 
    MyVector = c(beta)
    ####Split data
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0,];
    dat_obs = as.matrix(cbind("const" = rep(1, nrow(mydat)), mydat[,c("y4pred0")]));
    dat_obs[,1] = dat_obs[,1]/mydat$pred_obs_0; dat_obs[,2] = dat_obs[,2]/mydat$pred_obs_0
    tmp = rowSums(as.data.frame( t(t(dat_obs) * MyVector)) ) + mydat$h0;  tmp2 = sapply(tmp,function(x) (exp(-x) ));
    
    dat_obs = as.data.frame(cbind("const" = rep(1, nrow(mydat)), mydat[,c("y4pred0")]))
    dat_obs$term2 = tmp2
    dat_obs$const = dat_obs$const*dat_obs$term2;
    dat_obs$y4pred0 = dat_obs$y4pred0*dat_obs$term2;
    dat_obs$term2=NULL;
    
    U_obs = colSums(dat_obs)
    
    ####A3=0
    mydat = mytmpdat[mytmpdat$A_0!=mytmpdat$Q_0,]; 
    dat_miss = as.data.frame(cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred0")]))
    
    U_miss = -colSums(dat_miss)
    
    U = U_obs + U_miss
    c(U)
  }
  mybeta = c(0,0)
  #mybeta = c(0.005698283, -0.041822882, 0.017763672);
  #mybeta = c(-0.012722552, -0.008274138, 0.036790243)
  param0 = nleqslv(mybeta, model,control=list(allowSingular=FALSE))$x
  
  ### predict ###
  tmpdatweight =as.matrix(cbind("const" = rep(1, nrow(mytmpdat)), mytmpdat[,c("y4pred0")]));
  tmpdatweight[,1] = tmpdatweight[,1]/mytmpdat$pred_obs_0; tmpdatweight[,2] = tmpdatweight[,2]/mytmpdat$pred_obs_0
  tmp = rowSums(as.data.frame( t(t(tmpdatweight) * param0)) )  + mytmpdat$h0;  
  mytmpdat$weight0 = sapply(tmp,function(x) (plogis(x) ));
  
  ### repeat ###
  model=function(beta){ 
    MyVector = c(beta)
    ####Split data
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$Q_0==0 & mytmpdat$A_1==mytmpdat$Q_1,]; 
    dat_obs = as.matrix(cbind("const" = rep(1, nrow(mydat)), mydat[,c("y4pred1")]));
    dat_obs[,1] = dat_obs[,1]/(mydat$pred_obs_0*mydat$pred_obs_1); dat_obs[,2] = dat_obs[,2]/(mydat$pred_obs_0*mydat$pred_obs_1)
    tmp = rowSums(as.data.frame( t(t(dat_obs) * MyVector)) ) + mydat$h1;  tmp2 = sapply(tmp,function(x) (exp(-x) ));
    
    dat_obs = cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred1","weight0")])
    dat_obs$term2 = tmp2*(dat_obs$weight0)^-1
    dat_obs$const = dat_obs$const*dat_obs$term2; dat_obs$y4pred1 = dat_obs$y4pred1*dat_obs$term2;
    dat_obs$term2=NULL; dat_obs$weight0=NULL;
    
    U_obs = colSums(dat_obs)
    
    ####A3=0
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$Q_0==0 & mytmpdat$A_1!=mytmpdat$Q_1,]; 
    dat_miss = cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred1","weight0")])
    dat_miss$term2 = (dat_miss$weight0)^-1
    dat_miss$const = dat_miss$const*dat_miss$term2; dat_miss$y4pred1 = dat_miss$y4pred1*dat_miss$term2
    dat_miss$term2=NULL; dat_miss$weight0=NULL;    
    U_miss = -colSums(dat_miss)
    
    U = U_obs + U_miss
    c(U)
  }
  #mybeta = c(-0.03181454, -0.01592353, -0.05706790)
  #mybeta = c(-0.006886807, 0.002087208, 0.006178838)
  param1 = nleqslv(mybeta, model,control=list(allowSingular=FALSE))$x
  
  ### predict ###
  tmpdatweight =as.matrix(cbind("const" = rep(1, nrow(mytmpdat)), mytmpdat[,c("y4pred1")]));
  tmpdatweight[,1] = tmpdatweight[,1]/(mytmpdat$pred_obs_0*mytmpdat$pred_obs_1); tmpdatweight[,2] = tmpdatweight[,2]/(mytmpdat$pred_obs_0*mytmpdat$pred_obs_1)
  tmp = rowSums(as.data.frame( t(t(tmpdatweight) * param1)) ) + mytmpdat$h1;  
  mytmpdat$weight1 = sapply(tmp,function(x) (plogis(x) ));
  mytmpdat$weight1 = ifelse(mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_0==1, 1, mytmpdat$weight1)
  
  
  
  ### repeat ###
  model=function(beta){ 
    MyVector = c(beta)
    ####Split data
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_1==mytmpdat$Q_1 & mytmpdat$Q_1==0 & mytmpdat$A_2==mytmpdat$Q_2,]; 
    dat_obs = as.matrix(cbind("const" = rep(1, nrow(mydat)), mydat[,c("y4pred2")]));
    dat_obs[,1] = dat_obs[,1]/(mydat$pred_obs_0*mydat$pred_obs_1*mydat$pred_obs_2); dat_obs[,2] = dat_obs[,2]/(mydat$pred_obs_0*mydat$pred_obs_1*mydat$pred_obs_2)
    tmp = rowSums(as.data.frame( t(t(dat_obs) * MyVector)) ) + mydat$h2;  tmp2 = sapply(tmp,function(x) (exp(-x) ));
    
    dat_obs = cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred2","weight0","weight1")])
    dat_obs$term2 = tmp2*(dat_obs$weight0*dat_obs$weight1)^-1
    dat_obs$const = dat_obs$const*dat_obs$term2; dat_obs$y4pred2 = dat_obs$y4pred2*dat_obs$term2;
    dat_obs$term2=NULL; dat_obs$weight0=NULL; dat_obs$weight1=NULL;
    
    U_obs = colSums(dat_obs)
    
    ####A3=0
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_1==mytmpdat$Q_1 & mytmpdat$Q_1==0 & mytmpdat$A_2!=mytmpdat$Q_2,]; 
    dat_miss = cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred2","weight0","weight1")])
    dat_miss$term2 = (dat_miss$weight0*dat_miss$weight1)^-1
    dat_miss$const = dat_miss$const*dat_miss$term2; dat_miss$y4pred2 = dat_miss$y4pred2*dat_miss$term2
    dat_miss$term2=NULL; dat_miss$weight0=NULL; dat_miss$weight1=NULL;    
    U_miss = -colSums(dat_miss)
    
    U = U_obs + U_miss
    c(U)
  }
  mybeta = c(0,0)
  #mybeta = c(-0.198079591, 0.009709742, 0.358142519)
  #mybeta = c(0.021942593, -0.002201413, -0.034734934)
  param2 = nleqslv(mybeta, model,control=list(allowSingular=FALSE))$x
  
  ### predict ###
  tmpdatweight =as.matrix(cbind("const" = rep(1, nrow(mytmpdat)), mytmpdat[,c("y4pred2")]));
  tmpdatweight[,1] = tmpdatweight[,1]/(mytmpdat$pred_obs_0*mytmpdat$pred_obs_1*mytmpdat$pred_obs_2); tmpdatweight[,2] = tmpdatweight[,2]/(mytmpdat$pred_obs_0*mytmpdat$pred_obs_1*mytmpdat$pred_obs_2)
  tmp = rowSums(as.data.frame( t(t(tmpdatweight) * param2)) ) + mytmpdat$h2;  
  mytmpdat$weight2 = sapply(tmp,function(x) (plogis(x) ))
  mytmpdat$weight2 = ifelse(mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_1==mytmpdat$Q_1 & mytmpdat$A_1==1, 1, mytmpdat$weight2)
  
  
  ### repeat ###
  model=function(beta){ 
    MyVector = c(beta)
    ####Split data
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_1==mytmpdat$Q_1 & mytmpdat$A_2==mytmpdat$Q_2 & mytmpdat$Q_2==0 & mytmpdat$A_3==mytmpdat$Q_3,]; 
    dat_obs = as.matrix(cbind("const" = rep(1, nrow(mydat)), mydat[,c("y4pred3")]));
    dat_obs[,1] = dat_obs[,1]/(mydat$pred_obs_0*mydat$pred_obs_1*mydat$pred_obs_2*mydat$pred_obs_3); dat_obs[,2] = dat_obs[,2]/(mydat$pred_obs_0*mydat$pred_obs_1*mydat$pred_obs_2*mydat$pred_obs_3)
    tmp = rowSums(as.data.frame( t(t(dat_obs) * MyVector)) ) + mydat$h3;  tmp2 = sapply(tmp,function(x) (exp(-x) ));
    
    dat_obs = cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred3","weight0","weight1","weight2")])
    dat_obs$term2 = tmp2*(dat_obs$weight0*dat_obs$weight1*dat_obs$weight2)^-1
    dat_obs$const = dat_obs$const*dat_obs$term2; dat_obs$y4pred3 = dat_obs$y4pred3*dat_obs$term2;
    dat_obs$term2=NULL; dat_obs$weight0=NULL; dat_obs$weight1=NULL; dat_obs$weight2=NULL;
    U_obs = colSums(dat_obs)
    
    ####A3=0
    mydat = mytmpdat[mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_1==mytmpdat$Q_1 & mytmpdat$A_2==mytmpdat$Q_2 & mytmpdat$Q_2==0 & mytmpdat$A_3!=mytmpdat$Q_3,]; 
    dat_miss = cbind("const" = rep(1, nrow(mydat)),  mydat[,c("y4pred3","weight0","weight1","weight2")])
    dat_miss$term2 = (dat_miss$weight0*dat_miss$weight1*dat_miss$weight2)^-1
    dat_miss$const = dat_miss$const*dat_miss$term2; 
    dat_miss$y4pred3 = dat_miss$y4pred3*dat_miss$term2
    dat_miss$term2=NULL; dat_miss$weight0=NULL; dat_miss$weight1=NULL; dat_miss$weight2=NULL;     
    U_miss = -colSums(dat_miss)
    
    U = U_obs + U_miss
    c(U)
  }
  mybeta = c(0,0)
  #mybeta = c(0.39475759, -0.03394329, -0.49487539)
  #mybeta = c(0.399148956, -0.010877327, -0.513587440)
  param3 = nleqslv(mybeta, model,control=list(allowSingular=FALSE))$x
  
  
  ### predict ###
  tmpdatweight =as.matrix(cbind("const" = rep(1, nrow(mytmpdat)), mytmpdat[,c("y4pred3")]));
  tmpdatweight[,1] = tmpdatweight[,1]/(mytmpdat$pred_obs_0*mytmpdat$pred_obs_1*mytmpdat$pred_obs_2*mytmpdat$pred_obs_3); tmpdatweight[,2] = tmpdatweight[,2]/(mytmpdat$pred_obs_0*mytmpdat$pred_obs_1*mytmpdat$pred_obs_2*mytmpdat$pred_obs_3)
  tmp = rowSums(as.data.frame( t(t(tmpdatweight) * param3)) ) + mytmpdat$h3;  
  mytmpdat$weight3 = sapply(tmp,function(x) (plogis(x) ))
  mytmpdat$weight3 = ifelse(mytmpdat$A_0==mytmpdat$Q_0 & mytmpdat$A_1==mytmpdat$Q_1 & mytmpdat$A_2==mytmpdat$Q_2 & mytmpdat$A_2==1, 1, mytmpdat$weight3)
  
  tmpdata = mytmpdat
  
  tmpdata$pi0 = tmpdata$weight0
  tmpdata$pi1 = tmpdata$pi0*tmpdata$weight1
  tmpdata$pi2 = tmpdata$pi1*tmpdata$weight2
  tmpdata$pi3 = tmpdata$pi2*tmpdata$weight3
  
  tmp = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 &  
                  tmpdata$A_2==tmpdata$Q_2 &  tmpdata$A_3==tmpdata$Q_3,]; 
  tmpcount=0
  tmpcount = ifelse(any(tmp$pi0<0.00001,na.rm=T),1,tmpcount)
  tmpcount = ifelse(any(tmp$pi1<0.00001,na.rm=T),1,tmpcount)
  tmpcount = ifelse(any(tmp$pi2<0.00001,na.rm=T),1,tmpcount)
  tmpcount = ifelse(any(tmp$pi3<0.00001,na.rm=T),1,tmpcount)
  
  tmpdata$pi0 = ifelse(tmpdata$pi0<0.00001, 0.00001, tmpdata$pi0)
  tmpdata$pi1 = ifelse(tmpdata$pi1<0.00001, 0.00001, tmpdata$pi1)
  tmpdata$pi2 = ifelse(tmpdata$pi2<0.00001, 0.00001, tmpdata$pi2)
  tmpdata$pi3 = ifelse(tmpdata$pi3<0.00001, 0.00001, tmpdata$pi3)
  
  
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
  
  tmpcount2 = 0
  tmpcount2 = tmpcount2 + sum(is.na(tmp$meanytmp))
  
  myparam = c(meany4, tmpcount,tmpcount2)
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"WICE.csv")

stopCluster(cl)