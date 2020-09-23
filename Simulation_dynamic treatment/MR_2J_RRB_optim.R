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
  
  tmpdata$L0int = tmpdata$L1_0*tmpdata$L2_0
  tmpdata$L1int = tmpdata$L1_1*tmpdata$L2_1
  tmpdata$L2int = tmpdata$L1_2*tmpdata$L2_2
  tmpdata$L3int = tmpdata$L1_3*tmpdata$L2_3
  
  tmpdata$Q0int = tmpdata$Q_0*tmpdata$L2_0
  tmpdata$Q1int = tmpdata$Q_1*tmpdata$L2_1
  tmpdata$Q2int = tmpdata$Q_2*tmpdata$L2_2
  tmpdata$Q3int = tmpdata$Q_3*tmpdata$L2_3
  
  ##################
  ######time 3######
  ##################
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 &  
                    tmpdata$A_2==tmpdata$Q_2 &  tmpdata$A_3==tmpdata$Q_3,]; 
  myfit1 = glm(Y_4 ~ L1_3 + L2_3 + L1_3*L2_3 + Q_3 + L2_3*Q_3, family = binomial(), data = y4dat) ; 
  y4dat$y4pred1l = predict(myfit1, newdata = y4dat)
  mymodel = function(theta) {
    blah = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_2","L2_2","L2int","Q_2","Q2int")])) %*% theta;
    tmp1 = sum( (y4dat$pred_obs_3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred1l + blah)  ) );
    tmp2 = sum( (y4dat$pred_obs_3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L1_2 );
    tmp3 = sum( (y4dat$pred_obs_3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L2_2 );
    tmp4 = sum( (y4dat$pred_obs_3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L2int );
    tmp5 = sum( (y4dat$pred_obs_3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$Q_2 );
    tmp6 = sum( (y4dat$pred_obs_3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$Q2int );c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)}
  theta = c(0,0,0,0,0,0)
  myfit2  = nleqslv(theta, mymodel)$x
  y4dat$y4pred2l = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_2","L2_2","L2int","Q_2","Q2int")])) %*% myfit2 + y4dat$y4pred1l
  mymodel = function(theta) {
    blah = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% theta;
    tmp1 = sum( (y4dat$pred_obs_3*y4dat$pred_obs_2)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred2l + blah)  ) );
    tmp2 = sum( (y4dat$pred_obs_3*y4dat$pred_obs_2)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$L1_1 );
    tmp3 = sum( (y4dat$pred_obs_3*y4dat$pred_obs_2)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$L2_1 );
    tmp4 = sum( (y4dat$pred_obs_3*y4dat$pred_obs_2)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$L1int ); 
    tmp5 = sum( (y4dat$pred_obs_3*y4dat$pred_obs_2)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$Q_1 );
    tmp6 = sum( (y4dat$pred_obs_3*y4dat$pred_obs_2)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$Q1int );c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)}
  myfit3  = nleqslv(theta, mymodel)$x
  y4dat$y4pred3l = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% myfit3 + y4dat$y4pred2l;
  mymodel = function(theta) {
    blah = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_0","L2_0","L0int")])) %*% theta;
    tmp1 = sum( (y4dat$pi3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred3l + blah)  ) );
    tmp2 = sum( (y4dat$pi3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred3l + blah)  )*y4dat$L1_0 );
    tmp3 = sum( (y4dat$pi3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred3l + blah)  )*y4dat$L2_0 );
    tmp4 = sum( (y4dat$pi3)^-1*(y4dat$Y_4 - EXPIT(y4dat$y4pred3l + blah)  )*y4dat$L0int );c(tmp1,tmp2,tmp3,tmp4)}
  theta2 = c(0,0,0,0)
  myfit4  = nleqslv(theta2, mymodel)$x
  
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1 &  
                    tmpdata$A_2==tmpdata$Q_2,]; 
  y4dat$y4pred1 = predict(myfit1, newdata = y4dat, type="response"); 
  y4dat$y4pred1l = predict(myfit1, newdata = y4dat) 
  y4dat$y4pred2 = EXPIT(as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_2","L2_2","L2int","Q_2","Q2int")])) %*% myfit2 + y4dat$y4pred1l)
  y4dat$y4pred2l = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_2","L2_2","L2int","Q_2","Q2int")])) %*% myfit2 + y4dat$y4pred1l
  y4dat$y4pred3 = EXPIT(as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% myfit3 + y4dat$y4pred2l)
  y4dat$y4pred3l = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% myfit3 + y4dat$y4pred2l
  y4dat$y4pred4 = EXPIT(as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_0","L2_0","L0int")])) %*% myfit4 + y4dat$y4pred3l)
  y4dat$y4pred1 = ifelse(y4dat$Y_3==0,0,y4dat$y4pred1); y4dat$y4pred2 = ifelse(y4dat$Y_3==0,0,y4dat$y4pred2); 
  y4dat$y4pred3 = ifelse(y4dat$Y_3==0,0,y4dat$y4pred3); y4dat$y4pred4 = ifelse(y4dat$Y_3==0,0,y4dat$y4pred4);
  
  secondfit1 = glm(y4pred2 ~ L1_2 + L2_2 + L1_2*L2_2 + Q_2 + L2_2*Q_2, family = binomial(), data = y4dat) ; 
  y4dat$y4pred1l = predict(secondfit1, newdata = y4dat); 
  mymodel = function(theta) {
    blah = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% theta;
    tmp1 = sum( (y4dat$pred_obs_2)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  ) );
    tmp2 = sum( (y4dat$pred_obs_2)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L1_1 );
    tmp3 = sum( (y4dat$pred_obs_2)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L2_1 );
    tmp4 = sum( (y4dat$pred_obs_2)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L1int );
    tmp5 = sum( (y4dat$pred_obs_2)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$Q_1 );
    tmp6 = sum( (y4dat$pred_obs_2)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$Q1int ); c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)}
  seconfit2 = nleqslv(theta, mymodel)$x
  y4dat$y4pred2l = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% seconfit2 + y4dat$y4pred1l
  mymodel = function(theta) {
    blah = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_0","L2_0","L0int")])) %*% theta;
    tmp1 = sum( (y4dat$pi2)^-1*(y4dat$y4pred4 - EXPIT(y4dat$y4pred2l + blah)  ) );
    tmp2 = sum( (y4dat$pi2)^-1*(y4dat$y4pred4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$L1_0 );
    tmp3 = sum( (y4dat$pi2)^-1*(y4dat$y4pred4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$L2_0 );
    tmp4 = sum( (y4dat$pi2)^-1*(y4dat$y4pred4 - EXPIT(y4dat$y4pred2l + blah)  )*y4dat$L0int ); c(tmp1,tmp2,tmp3,tmp4)}
  seconfit3 = nleqslv(theta2, mymodel)$x
  
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0 & tmpdata$A_1==tmpdata$Q_1,]; 
  y4dat$y4pred1 = predict(secondfit1, newdata = y4dat, type="response");
  y4dat$y4pred1l = predict(secondfit1, newdata = y4dat)
  y4dat$y4pred2 = EXPIT(as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% seconfit2 + y4dat$y4pred1l)
  y4dat$y4pred2l = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_1","L2_1","L1int","Q_1","Q1int")])) %*% seconfit2 + y4dat$y4pred1l
  y4dat$y4pred3 = EXPIT(as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_0","L2_0","L0int")])) %*% seconfit3 + y4dat$y4pred2l) 
  y4dat$y4pred1 = ifelse(y4dat$Y_2==0,0,y4dat$y4pred1); y4dat$y4pred2 = ifelse(y4dat$Y_2==0,0,y4dat$y4pred2); y4dat$y4pred3 = ifelse(y4dat$Y_2==0,0,y4dat$y4pred3); 
  
  thirdfit1 = glm(y4pred2 ~ L1_1 + L2_1 + L1_1*L2_1 + Q_1 + L2_1*Q_1, family = binomial(), data = y4dat) ; 
  y4dat$y4pred1l = predict(thirdfit1, newdata = y4dat); 
  mymodel = function(theta) {
    blah = as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_0","L2_0","L0int")])) %*% theta;
    tmp1 = sum( (y4dat$pi1)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  ) );
    tmp2 = sum( (y4dat$pi1)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L1_0 );
    tmp3 = sum( (y4dat$pi1)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L2_0 );
    tmp4 = sum( (y4dat$pi1)^-1*(y4dat$y4pred3 - EXPIT(y4dat$y4pred1l + blah)  )*y4dat$L0int );c(tmp1,tmp2,tmp3,tmp4)}
  thirdfit2 = nleqslv(theta2, mymodel)$x
  
  y4dat = tmpdata[tmpdata$A_0==tmpdata$Q_0,]; 
  y4dat$y4pred1 = predict(thirdfit1, newdata = y4dat, type="response");
  y4dat$y4pred1l = predict(thirdfit1, newdata = y4dat)
  y4dat$y4pred2 = EXPIT(as.matrix(cbind("const" = rep(1, nrow(y4dat)), y4dat[,c("L1_0","L2_0","L0int")])) %*% thirdfit2 + y4dat$y4pred1l) 
  y4dat$y4pred1 = ifelse(y4dat$Y_1==0,0,y4dat$y4pred1);y4dat$y4pred2 = ifelse(y4dat$Y_1==0,0,y4dat$y4pred2); 
  
  finalfit = glm(y4pred2 ~ L1_0 + L2_0 + L1_0*L2_0, weight=1/pi0, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata; 
  y4dat$y4pred = predict(finalfit, newdata = y4dat, type="response"); 
  
  meany4tmp = c(mean(y4dat$y4pred))
  
  meany4 = (meany4tmp)
  #meany4
  
  myparam = cbind(meany4)
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"MR_2J_optim.csv")

stopCluster(cl)