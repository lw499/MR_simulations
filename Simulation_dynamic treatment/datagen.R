# Function to create survival data in long form
# It simulates ONE INDIVIDUAL's data at t = {0, 1, 2, ..., K-1} (outcomes realized at t = {1, 2, ..., K})
# Y is binary outcome
# Right censoring is present (C)
# A is binary treatment
# No competing risk
# 2 covariates (L1, L2)

# Inputs are i: id, df: data frame specified above, K: max number of time points
# Alphas, betas, thetas, etas, and sigma are user-provided parameters for data generating functions
datagen <- function(i, K, alpha0, alpha1, alpha2, alpha3, beta0, beta1, beta2,
                    theta0, theta1, theta2, theta3, theta4, kappa0, kappa1, kappa2, kappa3){
  id <- as.numeric(i)
  # Data at time 0
  L1 <- rbinom(1, 1, 0.5)
  L2 <- rbinom(1, 1, 0.5)
  A <- rbinom(1, 1, plogis(alpha0+alpha1*L2+alpha2*L1+alpha3*L1*L2))
  C <- 0
  Y <- rbinom(1, 1, plogis(theta0+theta1*A+theta2*L2+theta3*L1+theta4*L1*L2))
  #theta0+theta1*A0+theta2*L20+theta3*L10+theta4*L10*L20
  Q <- I(L1==0)
  temp <- data.frame(id = id, t0 = 0, L1, L2, A, C, Y, Q)
  
  if (temp$C==1){
    temp$Y <- NA
  }else if (temp$Y==1){
    for (j in 2:K){
      t0 <- j-1
      L1star <- as.numeric(rbinom(1, 1, plogis(kappa0+kappa1*temp$A[j-1]+kappa2*temp$L1[j-1]+kappa3*temp$L2[j-1])) )
      L2star <- as.numeric(rbinom(1, 1, plogis(beta0+beta1*temp$A[j-1]+beta2*L1star)) | temp$L2[j-1])
      #kappa0+kappa1*A0+kappa2*L10+kappa3*L20
      #(beta0+beta1*A0+beta2*L11) | L20)
      #temp_L1 <- c(temp$L1, L1star); temp_L2 <- c(temp$L2, L2star)
      Qstar = I(L1star==0) | temp$Q[j-1]
      Astar <- as.numeric(rbinom(1, 1, plogis(alpha0+alpha1*L2star+alpha2*L1star+alpha3*L1star*L2star)) | temp$A[j-1])
      #temp_A <- c(temp$A, Astar)
      Cstar <- 0
      if (Cstar==1){
        Ystar <- NA
        temp <- rbind(temp, c(id, t0, L1star, L2star, Astar, Cstar, Ystar, Qstar))
        break
      }
      else{
        Ystar <- rbinom(1, 1, plogis(theta0+theta1*Astar+theta2*L2star+theta3*L1star+theta4*L1star*L2star))
      }
      temp <- rbind(temp, c(id, t0, L1star, L2star, Astar, Cstar, Ystar, Qstar))
      if(Ystar!=1){break}
    }
  }
  return (temp)
}

#Qstar keeps track of the first time L1 becomes 0