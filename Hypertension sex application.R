# Continuous covariate

# derive var-cov matrix for each trial, then do a 2-stage IPDMA

# Lemonte paper, first order approximation to var-cov matrix for censored exp model
# - this is a 4 by 4 matrix for our model with an interaction term
# defined by X'WX = inv(I) - and we will calculate one for each trial

# Artur J. Lemonte (2020): Covariance matrix of maximum likelihood estimators
# in censored exponential regression models, Communications in Statistics - Theory and Methods,
# DOI: 10.1080/03610926.2020.1767142

# but the var-cov matrix is also inv(I*)/n where I* is the matrix based on the expected value

# As per logistic regression, we will generate a large dataset that matches 
# the characteristics of each trial 
# the values of each element for each patient, and then take their mean.
# This will give the expected values, which form the information matrix 
# And avoid having to specify a large matrix directly for X'WX

# run from START to END *
  
  ### START
DAT <- data.frame(n_C            = c(750,    199,    82,     2371,   3445,   1337,   2371,   131,    1139,   2297),
                  n_T            = c(780,    150,    90,     2427,   3546,   1314,   2365,   137,    1252,   2398),
                  age_mean_C     = c(42.36,  69.57,  74.11,  41.54,  45.17,  70.43,  71.54,  75.90,  66.77,  70.21),
                  age_sd_C       = c(5.34,   5.39,   8.69,   5.48,   5.86,   2.72,   6.68,   3.95,   5.67,   6.67),
                  age_mean_T     = c(42.17,  69.71,  72.64,  41.58,  45.38,  70.41,  71.64,  76.00,  66.42,  70.26),
                  age_sd_T       = c(5.39,   5.18,   7.99,   5.53,   6.00,   2.74,   6.72,   3.75,   5.34,   6.73),
                  percent_male_C = c(70.00,  37.19,  20.73,  53.73,  59.01,  41.81,  42.68,  24.43,  63.65,  33.83),
                  percent_male_T = c(69.36,  32.67,  25.56,  54.64,  58.83,  41.25,  43.72,  27.01,  65.02,  32.53),
                  events_C       = c(13,     28,     29,     82,     69,     199,    242,    7,      82,     137),
                  events_T       = c(9,      27,     32,     81,     73,     178,    213,    4,      61,     123),
                  trial          = c("ANBP", "COOP", "EWPH", "HDFP", "MRC1", "MRC2", "SHEP", "STOP", "SYCH", "SYSE"),
                  maxfup         = c(5.68,   10.34, 11.39,   5.00,   5.50,   8.39,   5.89,   5.18,   7.85,   8.11),
                  meanfup        = c(3.034,  4.580, 4.356,   4.915,  4.963,  5.770,  4.323,  1.955,  2.678,  2.493),
                  meanfup_T      = c(2.922,  4.426, 4.131,   4.911,  4.972,  5.749,  4.309,  1.926,  2.543,  2.486),
                  meanfup_C      = c(3.141,  4.785, 4.561,   4.920,  4.954,  5.791,  4.337,  1.983,  2.801,  2.500),
                  trialid        = c(1,      2,     3,       4,      5,      6,      7,      8,      9,      10))

DAT$rate_T <- DAT$events_T / (DAT$n_T * DAT$meanfup_T)
DAT$rate_C <- DAT$events_C / (DAT$n_C * DAT$meanfup_C)

#first define the  equation parameters for each study

# define vector alpha = log-rate in the control group
alpha <- log(DAT$rate_C)

# define vector of max follow up time
max_fup <- DAT$maxfup

# define vector of mean follow up time
mean_fup   <- DAT$meanfup
mean_fup_T <- DAT$meanfup_T
mean_fup_C <- DAT$meanfup_C

# define vector beta = overall treatment effect
beta <- log(DAT$rate_T / DAT$rate_C)

#* specify gamma = assumed prognostic effect of Z - assumed common for each trial here
gamma.vec <-  rep(log(1), nrow(DAT))
#gamma.vec <-  rep(log(1.5), nrow(DAT))
#gamma.vec <-  rep(log(1.25), nrow(DAT))

# specify assumed interaction - assumed common for each trial here
# say 1.3 for sex
# corresponds to lnHR of 0.262364

interaction <- rep(log(1.3), nrow(DAT))

# set up the sex variable
DAT$n_C_males <- DAT$percent_male_C * DAT$n_C / 100
DAT$n_C_males <- round(DAT$n_C_males, 0)

DAT$n_T_males <- DAT$percent_male_T * DAT$n_T / 100
DAT$n_T_males <- round(DAT$n_T_males, 0)

# total sample size of each trial
DAT$total = DAT$n_C + DAT$n_T

# total events of each trial

DAT$total_events <- DAT$events_C + DAT$events_T

# generate large dataset for the simulation
obs <- 1000000

# generate the number per group and covariate Z of interest; 
obs     <- 1000000
treat   <- matrix(0, ncol = nrow(DAT), nrow = obs)
control <- matrix(0, ncol = nrow(DAT), nrow = obs)
id      <- seq(1, obs, 1)
z_T     <- matrix(1,  ncol = nrow(DAT), nrow = obs)
z_C     <- matrix(1,  ncol = nrow(DAT), nrow = obs)
z       <- matrix(NA, ncol = nrow(DAT), nrow = obs)

for(i in 1:nrow(DAT)){
  x1 <- treat[,i]
  x1[id > obs * DAT$n_C[i] / DAT$total[i]] <- 1
  treat[,i] <- x1
  control[x1 == 0, i] <- 1

  x2 <- z_T[,i]
  x2[x1 == 1 & id > (obs * DAT$n_C[i] / DAT$total[i]) + (obs * DAT$n_T_males[i] / DAT$total[i] )] <- 0
  z_T[,i] <- x2
  z_C[x1 == 0 & id > obs * DAT$n_C_males[i] / DAT$total[i], i] <- 0

  z[,i] <- z_T[, i]
  z[x1 == 0, i] <- z_C[x1 == 0, i]
}

z_cent <- matrix(NA, ncol = nrow(DAT), nrow = obs)
for(i in 1:nrow(DAT)){
    z_cent[,i] <- scale(z[,i], scale = F)
}

# now generate the 4 by 4 matrix entries corresponding to ln(time) = alpha + beta*x + gamma*z + interaction*z*x
# for each study separately 

# So we want to derive X'WX, which is a 4 by 4 matrix 
# the elements of this are defined by x, z and z*x, and squared terms of these
# plus the W which is based on the probability of not being censored 


LP <- matrix(NA, ncol = nrow(DAT), nrow = obs)
W  <- matrix(NA, ncol = nrow(DAT), nrow = obs)
for(i in 1:nrow(DAT)){
	# Note that AFT models give same estimates as exponential reg model, but with opposite signs, and so I add negative value before each parameter
  LP[,i] <- -alpha[i] + (-beta[i]*treat[,i]) + (-gamma.vec[i] * z_cent[,i]) + (-interaction[i] * z_cent[,i] * treat[,i])

# Rather use original exp reg model framework
# gen LP`i' = alpha[1,`i'] + (beta[1,`i']*treat`i') + (gamma[1,`i'] * z`i'_cent) + (interaction[1,`i'] * z`i'_cent * treat`i')

# specifying W is tricky, and there are various options as follows - but mean does well on tests

# max f-up ... this will tend to give too low SEs
 # gen W`i' = 1 - exp(-max_fup[1,`i']*(exp(-LP`i')))
 
# based on overall prop of events - does ok I think, so may be an option
# corresponds to the type II censoring mentioned in the Lemonte paper
# gen W`i' = n_events[1,`i'] / n[1,`i']
 
# base on overall prop of events in each group - not much different to previous
# gen W`i' = events_T[1,`i'] / n_T[1,`i']
# replace  W`i' = events_C[1,`i'] / n_C[1,`i'] if control`i' == 1 
                           
                           # Use mean f-up - this makes more logical sense and tends to do well
                           # AFT approach
                           # gen W`i' = 1 - exp(-mean_fup[1,`i']*(exp(-LP`i')))

# exp reg model approach
# gen W`i' = mean_fup[1,`i']*(exp(LP`i'))

# AFT model - Kalbfleisch book approach - actually nearly equivalent to the Lemonte approach

W[,i] <- exp(log(mean_fup_T[i]) - (LP[,i]))
W[control[,i]==1,i] <- exp(log(mean_fup_C[i]) - (LP[control[,i]==1,i]))


M_11 <- matrix(ncol = nrow(DAT), nrow = obs)
M_12 <- matrix(ncol = nrow(DAT), nrow = obs)
M_13 <- matrix(ncol = nrow(DAT), nrow = obs)
M_14 <- matrix(ncol = nrow(DAT), nrow = obs)
M_21 <- matrix(ncol = nrow(DAT), nrow = obs)
M_22 <- matrix(ncol = nrow(DAT), nrow = obs)
M_23 <- matrix(ncol = nrow(DAT), nrow = obs)
M_24 <- matrix(ncol = nrow(DAT), nrow = obs)
M_31 <- matrix(ncol = nrow(DAT), nrow = obs)
M_32 <- matrix(ncol = nrow(DAT), nrow = obs)
M_33 <- matrix(ncol = nrow(DAT), nrow = obs)
M_34 <- matrix(ncol = nrow(DAT), nrow = obs)
M_41 <- matrix(ncol = nrow(DAT), nrow = obs)
M_42 <- matrix(ncol = nrow(DAT), nrow = obs)
M_43 <- matrix(ncol = nrow(DAT), nrow = obs)
M_44 <- matrix(ncol = nrow(DAT), nrow = obs)

  for(i in 1:nrow(DAT)){
    M_11[,i] <- W[,i]
    M_12[,i] <- treat[,i] * W[,i]
    M_13[,i] <- z_cent[,i] * W[,i]
    M_14[,i] <- treat[,i] * z_cent[,i] * W[,i]
    M_21[,i] <- M_12[,i]
    M_22[,i] <- treat[,i] * treat[,i] * W[,i]
    M_23[,i] <- treat[,i] * z_cent[,i] * W[,i]
    M_24[,i] <- treat[,i] * treat[,i] * z_cent[,i] * W[,i]
    M_31[,i] <- M_13[,i]
    M_32[,i] <- M_23[,i]
    M_33[,i] <- z_cent[,i] * z_cent[,i] * W[,i]
    M_34[,i] <- z_cent[,i] * z_cent[,i] * treat[,i] * W[,i]
    M_41[,i] <- M_14[,i]
    M_42[,i] <- M_24[,i]
    M_43[,i] <- M_34[,i]
    M_44[,i] <- z_cent[,i] * treat[,i] * z_cent[,i] * treat[,i] * W[,i]
  }
}

Imat <- array(dim = c(4, 4, nrow(DAT)))
for(i in 1:nrow(DAT)){
  Imat[1,1,] <- colMeans(M_11)
  Imat[1,2,] <- colMeans(M_12)
  Imat[1,3,] <- colMeans(M_13)
  Imat[1,4,] <- colMeans(M_14)
  Imat[2,1,] <- colMeans(M_21)
  Imat[2,2,] <- colMeans(M_22)
  Imat[2,3,] <- colMeans(M_23)
  Imat[2,4,] <- colMeans(M_24)
  Imat[3,1,] <- colMeans(M_31)
  Imat[3,2,] <- colMeans(M_32)
  Imat[3,3,] <- colMeans(M_33)
  Imat[3,4,] <- colMeans(M_34)
  Imat[4,1,] <- colMeans(M_41)
  Imat[4,2,] <- colMeans(M_42)
  Imat[4,3,] <- colMeans(M_43)
  Imat[4,4,] <- colMeans(M_44)
}



Imat_inv <- array(dim=c(4, 4, nrow(DAT)))
for(i in 1:nrow(DAT)){
  Imat_inv[,,i] <- solve(Imat[,,i])
}


# variance matrix of the parameters can be found by
# the 4,4 element is the asymptotic variance of the interaction 
vars <- vector(mode = 'numeric', length = nrow(DAT))
X    <- array(dim = c(4, 4, nrow(DAT)))
for(i in 1:nrow(DAT)){
  X[,,i]  <- Imat_inv[,,i] / DAT$total[i]
  vars[i] <- X[4,4,i]
}

se <-  sqrt(vars)
interact = interaction[i]

inv_var      <-  1/vars
summ_inv_var <- sum(inv_var)
ma_variance  <- 1/sum(inv_var)
ma_se        <- sqrt(1/sum(inv_var))

overall.var <- 1 / sum((c(vars))^-1)
overall.power <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[1] / sqrt(overall.var)) + pnorm(-qnorm(0.025, lower = F) - interaction[1] / sqrt(overall.var)))

power.ind <- vector(mode = 'numeric', length = nrow(DAT))
for(i in 1:10){
  power.ind[i] <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[i] / sqrt(vars[i])) + pnorm(-qnorm(0.025, lower = F) - interaction[i] / sqrt(vars[i])))
}

weight <- inv_var/summ_inv_var
percent_weight <- 100 * weight
percent_weight

data.frame(Study = 1:nrow(DAT), variance = vars, power = power.ind, weights = percent_weight)
overall.power

# Kalbfleisch book approach, prognostic effect of sex ln(1)
# interaction = .262364
# ma_variance = .01000047
# ma_se = .10000234
# lower = .06635942
# upper = .45836858
# power = .74652227
# power (%) = 74.652227


# Kalbfleisch book approach, prognostic effect of sex ln(1.5)
# interaction = .262364
# ma_variance = .00986575
# ma_se = .0993265
# lower = .06768406
# upper = .45704394
# power = .75220236
# power (%) = 75.220236


# Kab book approach, sex effect ln(1) and allowing for +1 and -1 to the mean followup for females and male
# interaction = .262364
# ma_variance = .01082942
# ma_se = .10406452
# lower = .05839754
# upper = .46633046
# power = .71266177
# power (%) = 71.266177


