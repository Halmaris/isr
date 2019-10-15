# Test Jarque–Bera for functional data

library(fda)
library(car) # Checking one-dimensional independence
library(MVN) # Multivariate normal tests (HZ, Mardia, Royston, energy)
library(normwhn.test) # Multivariate normal tests (DH)
require(moments) # Kurtosis & skewness
library(Hmisc) # To latex write
library(SuppDists) # For Johnson distribution
library(sn) # For skew normal distribution
library(gamlss.dist) # For the Sinh-Arcsinh (SHASH) distribution
# Jones & Pewsey (2009). Sinh-arcsinh distributions. Biometrika 96(4):761-780.

# Test Jarque–Bera (MJB_M & MJB_M*)
# Koizumi, Okamoto & Seo (2009). On Jarque–Bera tests for assessing multivariate normality. Journal of Statistics: Advances in Theory and Applications 1:207-220.
mjb.m = function(data){
  DNAME = deparse(substitute(data))
  data = scale(data, scale = FALSE) # Data centering
  data = as.matrix(data)
  n = dim(data)[1]
  p = dim(data)[2]
  S = ((n - 1) / n) * cov(data)
  D = data %*% solve(S) %*% t(data)
  bm1 = sum(D^3) / n^2
  bm2 = sum(diag((D^2))) / n
  zm1 = (n / 6) * bm1 * (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  zm2 = sqrt((n + 3) * (n + 5)) * ((n + 1) * bm2 - p * (p + 2) * (n - 1)) / 
    sqrt(8 * p * (p + 2) * (n - 3) * (n - p - 1) * (n - p + 1))
  df = p * (p + 1) * (p + 2) / 6 + 1
  stat1 = n * (bm1 / 6 + (bm2 - p * (p + 2))^2 / (8 * p * (p + 2)))
  stat2 = zm1 + zm2^2
  p.value1 = pchisq(stat1, df, lower.tail = FALSE)
  p.value2 = pchisq(stat2, df, lower.tail = FALSE)
  RVAL = list(statistic = c(MJB_M = stat1, MJB_M_star = stat2), 
              p.value = c(p.value1, p.value2), method = 'Multivariate Jarque-Bera Test (MJB_M, MJB_M*)', 
              data.name = DNAME)
  class(RVAL) = 'htest'
  return(RVAL)
}

sym.size = function(K = 75, N = 50){ # For random walk
  # K - length
  # N - number of observation per sample
  data = matrix(NA, nrow = N, ncol = K)
  for (i in 1:N)
    data[i,] = cumsum(rnorm(K))/sqrt(K)
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range)
  data.fd = smooth.basisPar(data.time, t(data), data.basis)$fd
  fpca.scores = pca.fd(data.fd, nharm = 4)$scores # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.single = pchisq(stat, 2 * p, lower.tail = FALSE)
  p.value.hz = mvn(fpca.scores, mvnTest = 'hz', desc = FALSE)$multivariateNormality$`p value`
  mardia.result = as.numeric(levels(mvn(fpca.scores, mvnTest = 'mardia', desc = FALSE)$multivariateNormality$`p value`))
  p.value.royston = mvn(fpca.scores, mvnTest = 'royston', desc = FALSE)$multivariateNormality$`p value`
  mjb.result = mjb.m(fpca.scores)
  p.value.energy = mvn(fpca.scores, mvnTest = 'energy', desc = FALSE)$multivariateNormality$`p value`
  dh1.result = mvn(fpca.scores, mvnTest = 'dh', desc = FALSE)$multivariateNormality$`p value`
  dh2.result = normality.test2(fpca.scores)
  return(list(p.value.single, p.value.hz, mardia.result[1], mardia.result[2], 
              p.value.royston, mjb.result$p.value[1], mjb.result$p.value[2], p.value.energy, dh1.result[1], 
              dh2.result[1]))
}

sym.power = function(K = 75, N = 50, df = 2) # For random walk
  # K - length
  # N - number of observation per sample
  # df - degree of freedom for t distribution
{
  data = matrix(NA, nrow = N, ncol = K)
  for (i in 1:N)
    data[i,] = cumsum(rt(K, df)) / sqrt(K)
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range)
  data.fd = smooth.basisPar(data.time, t(data), data.basis)$fd
  fpca.scores = pca.fd(data.fd, nharm = 4)$scores # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.single = pchisq(stat, 2 * p, lower.tail = FALSE)
  p.value.hz = mvn(fpca.scores, mvnTest = 'hz', desc = FALSE)$multivariateNormality$`p value`
  mardia.result = as.numeric(levels(mvn(fpca.scores, mvnTest = 'mardia', desc = FALSE)$multivariateNormality$`p value`))
  p.value.royston = mvn(fpca.scores, mvnTest = 'royston', desc = FALSE)$multivariateNormality$`p value`
  mjb.result = mjb.m(fpca.scores)
  p.value.energy = mvn(fpca.scores, mvnTest = 'energy', desc = FALSE)$multivariateNormality$`p value`
  dh1.result = mvn(fpca.scores, mvnTest = 'dh', desc = FALSE)$multivariateNormality$`p value`
  dh2.result = normality.test2(fpca.scores)
  return(list(p.value.single, p.value.hz, mardia.result[1], mardia.result[2], 
              p.value.royston, mjb.result$p.value[1], mjb.result$p.value[2], p.value.energy, dh1.result[1], 
              dh2.result[1]))
}

sym.size.reg = function(K = 75, N = 50, size = TRUE, df = 2) # For regression errors
  # K - length
  # N - number of observation per sample
  # size - if TRUE we calculate size, if FALSE we calculate power
{
  data.errors = matrix(NA, nrow = N, ncol = K) # Errors for Y(t)
  yt = matrix(NA, nrow = N, ncol = K)
  xt = matrix(NA, nrow = N, ncol = K)
  for (i in 1:N)
  {
    if (size) data.errors[i, ] = cumsum(rnorm(K)) / sqrt(K)
    else data.errors[i, ] = cumsum(rt(K, df)) / sqrt(K)
    yt[i, ] = (1:75) * i
    xt[i, ] = (1:75) * i
  }
  data = yt / 3 + data.errors # Division by 3 from calculation for psi(t, s) = ts and X(s) = ns.
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range)
  y.fd <<- smooth.basisPar(data.time, t(data), data.basis)$fd # Global export for fRegress
  x.fd <<- smooth.basisPar(data.time, t(xt), data.basis)$fd
  
  model.regression = fRegress(y.fd ~ x.fd) # Funcional regression model
  model.predict = predict(model.regression) # Prediction from functional regression model
  
  model.residuals = model.predict$fd - y.fd # Matrix of coefficients for residuals (functional object)
  
  fpca.scores = pca.fd(model.residuals, nharm = 4)$scores # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.single = pchisq(stat, 2 * p, lower.tail = FALSE)
  p.value.hz = mvn(fpca.scores, mvnTest = 'hz', desc = FALSE)$multivariateNormality$`p value`
  mardia.result = as.numeric(levels(mvn(fpca.scores, mvnTest = 'mardia', desc = FALSE)$multivariateNormality$`p value`))
  p.value.royston = mvn(fpca.scores, mvnTest = 'royston', desc = FALSE)$multivariateNormality$`p value`
  mjb.result = mjb.m(fpca.scores)
  p.value.energy = mvn(fpca.scores, mvnTest = 'energy', desc = FALSE)$multivariateNormality$`p value`
  dh1.result = mvn(fpca.scores, mvnTest = 'dh', desc = FALSE)$multivariateNormality$`p value`
  dh2.result = normality.test2(fpca.scores)
  rm(y.fd) #remove global obejects
  rm(x.fd)
  return(list(p.value.single, p.value.hz, mardia.result[1], mardia.result[2], 
              p.value.royston, mjb.result$p.value[1], mjb.result$p.value[2], p.value.energy, dh1.result[1], 
              dh2.result[1]))
}

parms <- JohnsonFit(c(0, 1, 0.4, 3), moment = 'quant') # We are looking for parametres of Johnson distribution
# Parameters of skew normal distribution
# alpha <- 10
# delta <- alpha / sqrt(1 + alpha^2)
# (skew <- (4 - pi) / 2 * (delta * sqrt(2 / pi))^3 / (1 - 2 * delta^2 / pi)^1.5)
# (ex.kurt <- 2 * (pi - 3) * (delta * sqrt(2 / pi))^3 / (1 - 2 * delta^2 / pi)^2)
# cat('Johnson distribution type:', parms$type)
# kurtosis >= skewness^2 + 1 (always)
gen.skew.kurt = function(K = 75, distr = 'normal'){
  # Function to generate time series with Z_t ~ distr; K - length
  Z = switch(distr,
              normal = c(0, 0, rnorm(K)),
              t4 = c(0, 0, rt(K, df = 4)),
              lognormal = c(0, 0, rlnorm(K)),
              chi2 = c(0, 0, rchisq(K, df = 2)),
              exponential = c(0, 0, rexp(K) - 1),
              johnson = c(0, 0, rJohnson(K, parms)), 
              sn = c(0, 0, rsn(n = K, xi = 0, omega = 1, alpha = 10)),
              shasho = c(0, 0, rSHASHo(K, mu = 0, sigma = 1, nu = .4, tau = 0.35)))
  X = numeric(K)
  for (i in (1:K))
    X[i] = Z[i] * i + 0.5 * Z[i + 1] * cos(pi * i) + 0.25 * Z[i + 2] * sin(pi * i)
  return(X)
}

# For IID errors
sym.iid.sk = function(K = 75, N = 50, distr = 'normal')
  # K - length
  # N - number of observation per sample
  # distr - alternative distribution
{
  data = matrix(NA, nrow = N, ncol = K)
  for (i in 1:N)
    data[i,] = gen.skew.kurt(K = K, distr = distr)
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range)
  data.fd = smooth.basisPar(data.time, t(data), data.basis)$fd
  fpca.scores = pca.fd(data.fd, nharm = 4)$scores # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.single = pchisq(stat, 2 * p, lower.tail = FALSE)
  p.value.hz = mvn(fpca.scores, mvnTest = 'hz', desc = FALSE)$multivariateNormality$`p value`
  mardia.result = as.numeric(levels(mvn(fpca.scores, mvnTest = 'mardia', desc = FALSE)$multivariateNormality$`p value`))
  p.value.royston = mvn(fpca.scores, mvnTest = 'royston', desc = FALSE)$multivariateNormality$`p value`
  mjb.result = mjb.m(fpca.scores)
  p.value.energy = mvn(fpca.scores, mvnTest = 'energy', desc = FALSE)$multivariateNormality$`p value`
  dh1.result = mvn(fpca.scores, mvnTest = 'dh', desc = FALSE)$multivariateNormality$`p value`
  dh2.result = normality.test2(fpca.scores)
  return(list(p.value.single, p.value.hz, mardia.result[1], mardia.result[2], 
              p.value.royston, mjb.result$p.value[1], mjb.result$p.value[2], p.value.energy, dh1.result[1],
              dh2.result[1]))
}

# For IID errors - final function to calculate size (for distr = 'normal') and power of tests
sym.results = function(K = 75, distr = 'normal')
{
  sym.50 = replicate(1000, sym.iid.sk(K = K, distr = distr))
  a1 = c(sum(sym.50[1,]<0.1), sum(sym.50[1,]<0.05), sum(sym.50[1,]<0.01)) # sizes of single test (JB)
  a2 = c(sum(sym.50[2,]<0.1), sum(sym.50[2,]<0.05), sum(sym.50[2,]<0.01)) # sizes of HZ test
  a3 = c(sum(sym.50[3,]<0.1), sum(sym.50[3,]<0.05), sum(sym.50[3,]<0.01)) # sizes of Mardia skew test
  a4 = c(sum(sym.50[4,]<0.1), sum(sym.50[4,]<0.05), sum(sym.50[4,]<0.01)) # sizes of Mardia kurt test
  a5 = c(sum(sym.50[5,]<0.1), sum(sym.50[5,]<0.05), sum(sym.50[5,]<0.01)) # sizes of Royston test
  a6 = c(sum(sym.50[6,]<0.1), sum(sym.50[6,]<0.05), sum(sym.50[6,]<0.01)) # sizes of MJB_M test
  a7 = c(sum(sym.50[7,]<0.1), sum(sym.50[7,]<0.05), sum(sym.50[7,]<0.01)) # sizes of MJB_M* test
  a8 = c(sum(sym.50[8,]<0.1), sum(sym.50[8,]<0.05), sum(sym.50[8,]<0.01)) # sizes of energy test
  a9 = c(sum(sym.50[9,]<0.1), sum(sym.50[9,]<0.05), sum(sym.50[9,]<0.01)) # sizes of DH1 test
  a10 = c(sum(sym.50[10,]<0.1), sum(sym.50[10,]<0.05), sum(sym.50[10,]<0.01)) # sizes of DH2 test
  result.sym.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
  sym.150 = replicate(1000, sym.iid.sk(K = K, distr = distr, N = 150))
  b1 = c(sum(sym.150[1,]<0.1), sum(sym.150[1,]<0.05), sum(sym.150[1,]<0.01)) # sizes of single test
  b2 = c(sum(sym.150[2,]<0.1), sum(sym.150[2,]<0.05), sum(sym.150[2,]<0.01)) # sizes of HZ test
  b3 = c(sum(sym.150[3,]<0.1), sum(sym.150[3,]<0.05), sum(sym.150[3,]<0.01)) # sizes of Mardia skew test
  b4 = c(sum(sym.150[4,]<0.1), sum(sym.150[4,]<0.05), sum(sym.150[4,]<0.01)) # sizes of Mardia kurt test
  b5 = c(sum(sym.150[5,]<0.1), sum(sym.150[5,]<0.05), sum(sym.150[5,]<0.01)) # sizes of Royston test
  b6 = c(sum(sym.150[6,]<0.1), sum(sym.150[6,]<0.05), sum(sym.150[6,]<0.01)) # sizes of MJB_M test
  b7 = c(sum(sym.150[7,]<0.1), sum(sym.150[7,]<0.05), sum(sym.150[7,]<0.01)) # sizes of MJB_M* test
  b8 = c(sum(sym.150[8,]<0.1), sum(sym.150[8,]<0.05), sum(sym.150[8,]<0.01)) # sizes of energy test
  b9 = c(sum(sym.150[9,]<0.1), sum(sym.150[9,]<0.05), sum(sym.150[9,]<0.01)) # sizes of DH1 test
  b10 = c(sum(sym.150[10,]<0.1), sum(sym.150[10,]<0.05), sum(sym.150[10,]<0.01)) # sizes of DH2 test
  result.sym.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
  sym.450 = replicate(1000, sym.iid.sk(K = K, distr = distr, N = 450))
  c1 = c(sum(sym.450[1,]<0.1), sum(sym.450[1,]<0.05), sum(sym.450[1,]<0.01)) # sizes of single test
  c2 = c(sum(sym.450[2,]<0.1), sum(sym.450[2,]<0.05), sum(sym.450[2,]<0.01)) # sizes of HZ test
  c3 = c(sum(sym.450[3,]<0.1), sum(sym.450[3,]<0.05), sum(sym.450[3,]<0.01)) # sizes of Mardia skew test
  c4 = c(sum(sym.450[4,]<0.1), sum(sym.450[4,]<0.05), sum(sym.450[4,]<0.01)) # sizes of Mardia kurt test
  c5 = c(sum(sym.450[5,]<0.1), sum(sym.450[5,]<0.05), sum(sym.450[5,]<0.01)) # sizes of Royston test
  c6 = c(sum(sym.450[6,]<0.1), sum(sym.450[6,]<0.05), sum(sym.450[6,]<0.01)) # sizes of MJB_M test
  c7 = c(sum(sym.450[7,]<0.1), sum(sym.450[7,]<0.05), sum(sym.450[7,]<0.01)) # sizes of MJB_M* test
  c8 = c(sum(sym.450[8,]<0.1), sum(sym.450[8,]<0.05), sum(sym.450[8,]<0.01)) # sizes of energy test
  c9 = c(sum(sym.450[9,]<0.1), sum(sym.450[9,]<0.05), sum(sym.450[9,]<0.01)) # sizes of DH1 test
  c10 = c(sum(sym.450[10,]<0.1), sum(sym.450[10,]<0.05), sum(sym.450[10,]<0.01)) # sizes of DH2 test
  result.sym.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
  result.sym.all = cbind(result.sym.50, result.sym.150, result.sym.450)
  row.names(result.sym.all) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
  return(result.sym.all)
}

# For regression errors
sym.size.power.regression = function(K = 75, N = 50, distr = 'normal')
  # K - length
  # N - number of observation per sample
  # distr - distribution of errors
{
  data.errors = matrix(NA, nrow = N, ncol = K) # Errors for Y(t)
  yt = matrix(NA, nrow = N, ncol = K)
  xt = matrix(NA, nrow = N, ncol = K)
  for (i in 1:N) {
    data.errors[i, ] = gen.skew.kurt(K = K, distr = distr)
    yt[i, ] = (1:75) * i
    xt[i, ] = (1:75) * i
  }
  
  data = yt / 3 + data.errors
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range)
  y.fd <<- smooth.basisPar(data.time, t(data), data.basis)$fd # Global export for fRegress
  x.fd <<- smooth.basisPar(data.time, t(xt), data.basis)$fd
  
  model.regression = fRegress(y.fd ~ x.fd) # Funcional regression model
  model.predict = predict(model.regression) # Prediction from functional regression model
  
  model.residuals = model.predict$fd - y.fd # Matrix of coefficients for residuals (functional object)
  
  fpca.scores = pca.fd(model.residuals, nharm = 4)$scores # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.single = pchisq(stat, 2 * p, lower.tail = FALSE)
  p.value.hz = mvn(fpca.scores, mvnTest = 'hz', desc = FALSE)$multivariateNormality$`p value`
  mardia.result = as.numeric(levels(mvn(fpca.scores, mvnTest = 'mardia', desc = FALSE)$multivariateNormality$`p value`))
  p.value.royston = mvn(fpca.scores, mvnTest = 'royston', desc = FALSE)$multivariateNormality$`p value`
  mjb.result = mjb.m(fpca.scores)
  p.value.energy = mvn(fpca.scores, mvnTest = 'energy', desc = FALSE)$multivariateNormality$`p value`
  dh1.result = mvn(fpca.scores, mvnTest = 'dh', desc = FALSE)$multivariateNormality$`p value`
  dh2.result = normality.test2(fpca.scores)
  rm(y.fd) #remove global obejects
  rm(x.fd)
  return(list(p.value.single, p.value.hz, mardia.result[1], mardia.result[2], 
              p.value.royston, mjb.result$p.value[1], mjb.result$p.value[2], p.value.energy, dh1.result[1], 
              dh2.result[1]))
}

# For regression errors - final function to calculate size (for distr = 'normal') and power of tests
sym.results.regression = function(K = 75, distr = 'normal')
{
  sym.50 = replicate(1000, sym.size.power.regression(K = K, distr = distr))
  a1 = c(sum(sym.50[1,]<0.1), sum(sym.50[1,]<0.05), sum(sym.50[1,]<0.01)) # sizes of single test (JB)
  a2 = c(sum(sym.50[2,]<0.1), sum(sym.50[2,]<0.05), sum(sym.50[2,]<0.01)) # sizes of HZ test
  a3 = c(sum(sym.50[3,]<0.1), sum(sym.50[3,]<0.05), sum(sym.50[3,]<0.01)) # sizes of Mardia skew test
  a4 = c(sum(sym.50[4,]<0.1), sum(sym.50[4,]<0.05), sum(sym.50[4,]<0.01)) # sizes of Mardia kurt test
  a5 = c(sum(sym.50[5,]<0.1), sum(sym.50[5,]<0.05), sum(sym.50[5,]<0.01)) # sizes of Royston test
  a6 = c(sum(sym.50[6,]<0.1), sum(sym.50[6,]<0.05), sum(sym.50[6,]<0.01)) # sizes of MJB_M test
  a7 = c(sum(sym.50[7,]<0.1), sum(sym.50[7,]<0.05), sum(sym.50[7,]<0.01)) # sizes of MJB_M* test
  a8 = c(sum(sym.50[8,]<0.1), sum(sym.50[8,]<0.05), sum(sym.50[8,]<0.01)) # sizes of energy test
  a9 = c(sum(sym.50[9,]<0.1), sum(sym.50[9,]<0.05), sum(sym.50[9,]<0.01)) # sizes of DH1 test
  a10 = c(sum(sym.50[10,]<0.1), sum(sym.50[10,]<0.05), sum(sym.50[10,]<0.01)) # sizes of DH2 test
  result.sym.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
  sym.150 = replicate(1000, sym.size.power.regression(K = K, distr = distr, N = 150))
  b1 = c(sum(sym.150[1,]<0.1), sum(sym.150[1,]<0.05), sum(sym.150[1,]<0.01)) # sizes of single test
  b2 = c(sum(sym.150[2,]<0.1), sum(sym.150[2,]<0.05), sum(sym.150[2,]<0.01)) # sizes of HZ test
  b3 = c(sum(sym.150[3,]<0.1), sum(sym.150[3,]<0.05), sum(sym.150[3,]<0.01)) # sizes of Mardia skew test
  b4 = c(sum(sym.150[4,]<0.1), sum(sym.150[4,]<0.05), sum(sym.150[4,]<0.01)) # sizes of Mardia kurt test
  b5 = c(sum(sym.150[5,]<0.1), sum(sym.150[5,]<0.05), sum(sym.150[5,]<0.01)) # sizes of Royston test
  b6 = c(sum(sym.150[6,]<0.1), sum(sym.150[6,]<0.05), sum(sym.150[6,]<0.01)) # sizes of MJB_M test
  b7 = c(sum(sym.150[7,]<0.1), sum(sym.150[7,]<0.05), sum(sym.150[7,]<0.01)) # sizes of MJB_M* test
  b8 = c(sum(sym.150[8,]<0.1), sum(sym.150[8,]<0.05), sum(sym.150[8,]<0.01)) # sizes of energy test
  b9 = c(sum(sym.150[9,]<0.1), sum(sym.150[9,]<0.05), sum(sym.150[9,]<0.01)) # sizes of DH1 test
  b10 = c(sum(sym.150[10,]<0.1), sum(sym.150[10,]<0.05), sum(sym.150[10,]<0.01)) # sizes of DH2 test
  result.sym.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
  sym.450 = replicate(1000, sym.size.power.regression(K = K, distr = distr, N = 450))
  c1 = c(sum(sym.450[1,]<0.1), sum(sym.450[1,]<0.05), sum(sym.450[1,]<0.01)) # sizes of single test
  c2 = c(sum(sym.450[2,]<0.1), sum(sym.450[2,]<0.05), sum(sym.450[2,]<0.01)) # sizes of HZ test
  c3 = c(sum(sym.450[3,]<0.1), sum(sym.450[3,]<0.05), sum(sym.450[3,]<0.01)) # sizes of Mardia skew test
  c4 = c(sum(sym.450[4,]<0.1), sum(sym.450[4,]<0.05), sum(sym.450[4,]<0.01)) # sizes of Mardia kurt test
  c5 = c(sum(sym.450[5,]<0.1), sum(sym.450[5,]<0.05), sum(sym.450[5,]<0.01)) # sizes of Royston test
  c6 = c(sum(sym.450[6,]<0.1), sum(sym.450[6,]<0.05), sum(sym.450[6,]<0.01)) # sizes of MJB_M test
  c7 = c(sum(sym.450[7,]<0.1), sum(sym.450[7,]<0.05), sum(sym.450[7,]<0.01)) # sizes of MJB_M* test
  c8 = c(sum(sym.450[8,]<0.1), sum(sym.450[8,]<0.05), sum(sym.450[8,]<0.01)) # sizes of energy test
  c9 = c(sum(sym.450[9,]<0.1), sum(sym.450[9,]<0.05), sum(sym.450[9,]<0.01)) # sizes of DH1 test
  c10 = c(sum(sym.450[10,]<0.1), sum(sym.450[10,]<0.05), sum(sym.450[10,]<0.01)) # sizes of DH2 test
  result.sym.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
  result.sym.all = cbind(result.sym.50, result.sym.150, result.sym.450)
  row.names(result.sym.all) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
  return(result.sym.all)
}