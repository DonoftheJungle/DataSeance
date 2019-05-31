library(MASS)
library(glmnet)
library(dmutate)
library(knockoff)
library(ggplot2)
library(reshape2)

set.seed(1)

siglevls = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 8, 10)

trials = 100

p = 10000
p_0 = 50
k = 1000
n = 400

#W1 = matrix(runif(10000, 1, 10), p, k)
#W2 = matrix(runif(10000, -10, 10), p, k)
#S1 = W1 %*% t(W1) + diag(runif(p, 1, 10))
#S2 = W2 %*% t(W2) + diag(runif(p, 1, 10))
#S = S1 + 2*S2

#C = cov2cor(S)
#C2 = cov2cor(S2)
#C3 = diag(1, p)

mu = rep(0, p)

#Generate design matrix and normalize columns to unit euclidean length
#X = mvrnorm(n, mu, S)
#X = rlmvnorm(n, mu, C3)
L2nrms = diag(t(X)%*% X)
X = X %*% (diag(1/sqrt(L2nrms)))

#generate true variable set
beta = rep(0, p)
sigidx = sample(500, 50)


#function to compute sign-restricted knockoff statistics
k_stat = function(X, X_k, y) stat.glmnet_coefdiff(X, X_k, y, upper.limits = upperlim, lower.limits=lowerlim)

#function to generate knockoffs with either data-splitting or data-recycling
knockoffs_split = function(X) create.fixed(X, method='sdp')

knockoffs_recycle = function(X) {
  Xk0 = create.fixed(X, method='sdp')
  Xkrec = rbind(Xscreen, Xk0$Xk)
  Xrec = rbind(Xscreen, X)
  y = NULL
  rec_knocks = list(Xrec, Xkrec, y)
  attr(rec_knocks, "names") <- c("X", "Xk", "y")
  attr(rec_knocks, "class") <- "knockoff.variables"
  return(rec_knocks)
}

siglevls_results = matrix(0, ncol = 6, nrow=length(siglevls))

for (j in 1:length(siglevls)){

sigmag = siglevls[j]
beta[sigidx] = sigmag
#randomly distribute signs of true coefficients
signs = 2*rbinom(p, 1, 0.5) - 1
beta = beta*signs


dirFDRs = rep(0, trials)
beta_power = rep(0, trials)
dirbeta_power = rep(0, trials)
restr_beta_power = rep(0, trials)
restr_dirbeta_power = rep(0, trials)

for (k in 1:trials){
  
  #Generate Response Data
  Y = X %*% beta + rnorm(n)

  #split data to perform initial variable screening
  n_0 = (3/8)*n
  n_1 = n - n_0
  screenidx = sample(n, n_0)

  Yscreen = Y[screenidx]
  Yfilter = Y[-screenidx]

  Xscreen = X[screenidx, ]
  Xfilter = X[-screenidx, ]

  #perform variable screening using LASSO

  scrn = glmnet(Xscreen, Yscreen, family = 'gaussian', alpha=1, lambda.min.ratio = 0.001)
  #Of all models with fewer than n_1/2 variables, choose the largest
  modelidx = max(which(colSums(scrn$beta != 0) < (n_1/2)))
  beta_s0 = scrn$beta[, modelidx]
  s0 = as.numeric(which(beta_s0 != 0))

  #select screened variables and get initial sign estimates
  Xfilter = Xfilter[, s0]
  Xscreen = Xscreen[, s0]
  beta_s0 = as.numeric(beta_s0[s0])
  sign_hat = sign(beta_s0)

  #Create bounds for knockoff estimates using initial sign estimates from screening
  upperlim = rep(Inf, length(s0))
  upperlim = replace(upperlim, sign_hat<0, 0)
  upperlim = c(upperlim, upperlim)
  lowerlim = rep(-Inf, length(s0))
  lowerlim = replace(lowerlim, sign_hat>0, 0)
  lowerlim = c(lowerlim, lowerlim)
  
  #run knockoff filter on reduced model
  results = knockoff.filter(Xfilter, Yfilter, knockoffs=knockoffs_split, statistic=k_stat, fdr = 0.1)
  koidx = results$selected
  sign_ko = sign_hat[koidx]
  
  #compute a partial regression on the original data, but only using the variables from the screened set
  #will be used to compute dirFDR
  lmpart = lm(Yfilter ~ Xfilter)
  sign_part = as.numeric(sign(lmpart$coefficients))
  #remove the intercept term
  sign_part = sign_part[2:length(sign_part)]
  #extract the signs of the 'true' partial coefficients corresponding to the variables chosen by knockof filter
  sign_part_true = sign_part[koidx]
  
  #Count the number of directional errors among the coefficients from the reduced model s0, and compute dirFDR
  err_dir = sum(sign_part_true != sign_ko)
  dirFDRs[k] = err_dir/max(length(koidx), 1)
  
  #Retrieve the indices of the final set of selected variables
  select_idx = s0[koidx]
  #compute the power and directional power relative to the full parameter vector
  true_disc = intersect(select_idx, sigidx)
  beta_power[k] = length(true_disc)/p_0
  dirtrue_disc = sum(sign(beta[select_idx]) == sign_ko)
  dirbeta_power[k] = dirtrue_disc/p_0
  
  #compute the power and directional power restricted to the screened model
  restr_beta_power[k] = length(true_disc)/max(length(intersect(s0, sigidx)), 1)
  restr_dirbeta_power[k] = dirtrue_disc/max(length(intersect(s0, sigidx)), 1)
}

avg_dirFDR = mean(dirFDRs)
avg_beta_power = mean(beta_power)
avg_dirbeta_power = mean(dirbeta_power)
avg_restr_beta_power = mean(restr_beta_power)
avg_restr_dirbeta_power = mean(restr_dirbeta_power)
avg_results = c(sigmag, avg_dirFDR, avg_beta_power, avg_dirbeta_power, avg_restr_beta_power, avg_restr_dirbeta_power)
print(sprintf("Average dirFDR: %.3f  Average Power: %.3f   Average Directional Power: %.3f", 
              avg_dirFDR, avg_beta_power, avg_dirbeta_power))
print(sprintf("Average Restricted Power: %.3f   Average Restricted Directional Power: %.3f",
              avg_restr_beta_power, avg_restr_dirbeta_power))

siglevls_results[j, ] = avg_results
}

results_frame = as.data.frame(siglevls_results)
names(results_frame) = c("Signal", "dirFDR", "Power", "dirPower", "Restricted Power", "Restricted dirPower")
meltframe = melt(results_frame, id="Signal")
pl = ggplot(meltframe, aes(x = Signal, y = value, group = variable, colour = variable)) + 
    geom_line(aes(linetype=variable))
pl + scale_color_manual(values = rainbow(5))