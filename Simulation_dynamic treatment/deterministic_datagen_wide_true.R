set.seed(5678)

#L2: Indicator of AIDS
#L1: Indicator of high CD4 count 
#A: treatment

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
N <- 10000000
Ndat <- 4
beta0=-2; beta1=-2; beta2=-1;
##If last treatment is 1, lower chances AIDS (beta1 is negative)
#If current CD4 value is high, then chances of AIDS is low (beta2 negative)
theta0=-1; theta1=2; theta2=-2; theta3=2; theta4=2;
#If last treatment is 1, higher chances survival (theta1 is positive)
#If current AIDS, probability of survival lower (theta2 negative)
#If current CD4 count high, probability of survival higher (theta3 positive)
kappa0 = 0; kappa1 = 1; kappa2 = 1; kappa3 = -1; 
# If last treatment is 1, then higher chances CD4 high (kappa1 positive)
# If last L1 is 1 then there's a higher chance it will be high (kappa2 is positive)
# If last AIDS status is 1, then lower chance that L1 is high (kappa3 is negative)

L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Y5 <- L15 <- L25 <- A5 <- Y6 <- L16 <- L26 <- A6 <- Y7 <- as.numeric(rep(NA, N))
ids <- as.list(1:N)
L10 <- rbinom(N, 1, 0.5)
L20 <- rbinom(N, 1, 0.5)
#A0 <- rep(1, N)
A0 = I(L10==0)
Y1 <- plogis(theta0+theta1*A0+theta2*L20+theta3*L10+theta4*L10*L20)

L11 = as.numeric(rexpit(kappa0+kappa1*A0+kappa2*L10+kappa3*L20))
L21 = as.numeric(rexpit(beta0+beta1*A0+beta2*L11) | L20)
#A1 = rep(1,N)
A1 = I(L11==0) | A0 # If L11 is 0 (threshold is not crossed), then A1 is 0
Y2 = plogis((theta0+theta1*A1+theta2*L21+theta3*L11+theta4*L11*L21))

L12 = as.numeric(rexpit(kappa0+kappa1*A1+kappa2*L11+kappa3*L21))
L22 = as.numeric(rexpit(beta0+beta1*A1+beta2*L12) | L21)
#A2 = rep(1, N)
A2 = I(L12==0) | A1
Y3 = plogis((theta0+theta1*A2+theta2*L22+theta3*L12+theta4*L12*L22))

L13 = as.numeric(rexpit(kappa0+kappa1*A2+kappa2*L12+kappa3*L22))
L23 = as.numeric(rexpit(beta0+beta1*A2+beta2*L13) | L22)
#A3 = rep(1, N)
A3 = I(L13==0) | A2
Y4 = plogis((theta0+theta1*A3+theta2*L23+theta3*L13+theta4*L13*L23))

surv1 = mean(Y1)
surv2 = mean((Y1)*(Y2))
surv3 = mean((Y1)*(Y2)*(Y3))
surv4 = mean((Y1)*(Y2)*(Y3)*(Y4))

c(surv1, surv2, surv3, surv4)
#static: 
#dynamic (L1): 0.6155334 0.4261114 0.3258812 0.2661737  0.2263013  #start ART when CD4 first falls below low
