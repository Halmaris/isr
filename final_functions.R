# Test Jarque–Bera for functional data

library(fda) # For funcional data methods
library(MVN) # Multivariate normal tests (HZ, Mardia, Royston)
library(normwhn.test) # Multivariate normal tests (DH)
library(moments) # Kurtosis & skewness
library(Hmisc) # To latex write
library(sn) # For skew normal distribution

# Test Jarque–Bera (MJB_M)
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
  # zm1 = (n / 6) * bm1 * (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  df = p * (p + 1) * (p + 2) / 6 + 1
  stat1 = n * (bm1 / 6 + (bm2 - p * (p + 2))^2 / (8 * p * (p + 2)))
  p.value1 = pchisq(stat1, df, lower.tail = FALSE)
  RVAL = list(statistic = stat1, 
              p.value = p.value1, method = 'Multivariate Jarque-Bera Test (MJB_M)', 
              data.name = DNAME)
  class(RVAL) = 'htest'
  return(RVAL)
}

# Function to generate data
gen.skew.kurt = function(K = 75, distr = 'normal'){
  # Function to generate time series with Z_t ~ distr; K - length
  Z = switch(distr,
              normal = c(0, 0, rnorm(K)),
              t4 = c(0, 0, rt(K, df = 4)),
              exponential = c(0, 0, rexp(K) - 1),
              sn = c(0, 0, rsn(n = K, xi = 0, omega = 1, alpha = 10)))
  X = numeric(K)
  for (i in (1:K))
    X[i] = Z[i] * i + 0.5 * Z[i + 1] * cos(pi * i) + 0.25 * Z[i + 2] * sin(pi * i)
  return(X)
}

# For IID errors
sym.iid.sk = function(K = 75, N = 50, distr = 'normal', percent = 0.85, param = TRUE)
  # K - length
  # N - number of observation per sample
  # distr - alternative distribution
  # percent - percent of variation to explain
  # param - T/F; TRUE means 5.1 process; FALSE means 5.2 process
{
  data = matrix(NA, nrow = N, ncol = K)
  if (param) {
    for (i in 1:N) {
      Z = switch(distr,
                  normal = rnorm(K),
                  t4 = rt(K, df = 4),
                  exponential = rexp(K) - 1,
                  sn = rsn(n = K, xi = 0, omega = 1, alpha = 10))
      data[i, ] = cumsum(Z) / sqrt(K)  
    }
  } 
  else {
    for (i in 1:N)
      data[i, ] = gen.skew.kurt(K = K, distr = distr)
  }
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range, nbasis = 10)
  data.fd = smooth.basisPar(data.time, t(data), data.basis)$fd
  fpca.model = pca.fd(data.fd, nharm = 10) # FPCA model
  count <- min(which(cumsum(fpca.model$varprop) > percent)) # Dimension of functional data projection
  count <- ifelse(count == 1, 2, count)
  fpca.scores <- fpca.model$scores[, 1:count] # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.single = pchisq(stat, 2 * p, lower.tail = FALSE)
  mjb.result = mjb.m(fpca.scores)
  dh1.result = normality.test1(fpca.scores)
  dh2.result = normality.test2(fpca.scores)
  return(list(p.value.single, mjb.result$p.value, dh1.result[1], dh2.result[1], p = count))
}

# For IID errors - final function to calculate size (for distr = 'normal') and power of tests
sym.results = function(K = 75, distr = 'normal', counter = 1000, percent = 0.85, param = TRUE)
  # counter - number of replications
  # percent - percent of variation to explain
  # param - T/F; TRUE means 5.1 process; FALSE means 5.2 process
{
  sym.50 = replicate(counter, sym.iid.sk(K = K, distr = distr, percent = percent, param = param))
  a1 = c(sum(sym.50[1, ] < 0.1), sum(sym.50[1, ] < 0.05), sum(sym.50[1, ] < 0.01)) # JB test
  a2 = c(sum(sym.50[2, ] < 0.1), sum(sym.50[2, ] < 0.05), sum(sym.50[2, ] < 0.01)) # MJB_M test
  a3 = c(sum(sym.50[3, ] < 0.1), sum(sym.50[3, ] < 0.05), sum(sym.50[3, ] < 0.01)) # DH1 test
  a4 = c(sum(sym.50[4, ] < 0.1), sum(sym.50[4, ] < 0.05), sum(sym.50[4, ] < 0.01)) # DH2 test
  result.sym.50 = rbind(a1, a2, a3, a4)
  sym.150 = replicate(counter, sym.iid.sk(K = K, distr = distr, N = 150, percent = percent, param = param))
  b1 = c(sum(sym.150[1, ] < 0.1), sum(sym.150[1, ] < 0.05), sum(sym.150[1, ] < 0.01)) # JB test
  b2 = c(sum(sym.150[2,] < 0.1), sum(sym.150[2, ] < 0.05), sum(sym.150[2, ] < 0.01)) # MJB_M test
  b3 = c(sum(sym.150[3,] < 0.1), sum(sym.150[3, ] < 0.05), sum(sym.150[3, ] < 0.01)) # DH1 test
  b4 = c(sum(sym.150[4,] < 0.1), sum(sym.150[4, ] < 0.05), sum(sym.150[4, ] < 0.01)) # DH2 test
  result.sym.150 = rbind(b1, b2, b3, b4)
  sym.450 = replicate(counter, sym.iid.sk(K = K, distr = distr, N = 450, percent = percent, param = param))
  c1 = c(sum(sym.450[1, ] < 0.1), sum(sym.450[1, ] < 0.05), sum(sym.450[1, ] < 0.01)) # JB test
  c2 = c(sum(sym.450[2, ] < 0.1), sum(sym.450[2, ] < 0.05), sum(sym.450[2, ] < 0.01)) # MJB_M test
  c3 = c(sum(sym.450[3, ] < 0.1), sum(sym.450[3, ] < 0.05), sum(sym.450[3, ] < 0.01)) # DH1 test
  c4 = c(sum(sym.450[4, ] < 0.1), sum(sym.450[4, ] < 0.05), sum(sym.450[4, ] < 0.01)) # DH2 test
  result.sym.450 = rbind(c1, c2, c3, c4)
  result.sym.all = cbind(result.sym.50, result.sym.150, result.sym.450)
  row.names(result.sym.all) = c('Jarque-Bera',  'MJB_M', 'Doornik-Hansen', 'Lobato-Velasco')
  standard.error <- sqrt(result.sym.all / counter * (1 - result.sym.all / counter) / counter)
  p <- list(unlist(sym.50[5, ]), unlist(sym.150[5, ]), unlist(sym.450[5, ]))
  return(list(round(result.sym.all / counter, 3), round(standard.error, 3), p))
}

# For regression errors
sym.regression.sk = function(K = 75, N = 50, distr = 'normal', percent = 0.85, param = TRUE)
  # K - length
  # N - number of observation per sample
  # distr - alternative distribution
  # percent - percent of variation to explain
  # param - T/F; TRUE means 5.1 process; FALSE means 5.2 process
{
  data.errors = matrix(NA, nrow = N, ncol = K) # Errors for Y(t)
  yt = matrix(NA, nrow = N, ncol = K)
  xt = matrix(NA, nrow = N, ncol = K)
  
  if (param) {
    for (i in 1:N) {
      Z = switch(distr,
                  normal = rnorm(K),
                  t4 = rt(K, df = 4),
                  exponential = rexp(K) - 1,
                  sn = rsn(n = K, xi = 0, omega = 1, alpha = 10))
      data.errors[i, ] = cumsum(Z) / sqrt(K)  
      yt[i, ] = (1:K) * i
      xt[i, ] = (1:K) * i
    }
  } 
  else {
    for (i in 1:N) {
      data.errors[i, ] = gen.skew.kurt(K = K, distr = distr)
      yt[i, ] = (1:K) * i
      xt[i, ] = (1:K) * i
    }
  }
  
  data = yt / (3 * K) + data.errors
  data.time = seq(0, 1, length = K)
  data.range = c(0, 1)
  data.basis = create.bspline.basis(data.range, nbasis = 4)
  y.fd <<- smooth.basisPar(data.time, t(data), data.basis)$fd # Global export for fRegress
  x.fd <<- smooth.basisPar(data.time, t(xt), data.basis)$fd
  
  model.regression = fRegress(y.fd ~ x.fd) # Funcional regression model
  model.predict = predict(model.regression) # Prediction from functional regression model
  
  model.residuals = model.predict$fd - y.fd # Matrix of coefficients for residuals (functional object)
  
  fpca.model = pca.fd(model.residuals, nharm = 4) # FPCA model
  count <- min(which(cumsum(fpca.model$varprop) > percent)) # Dimension of functional data projection
  count <- ifelse(count == 1, 2, count)
  fpca.scores <- fpca.model$scores[, 1:count] # FPCA scores
  data.kurt = sapply(data.frame(fpca.scores), kurtosis)
  data.skewness = sapply(data.frame(fpca.scores), skewness)
  stat = N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p = length(data.skewness)
  p.value.jb = pchisq(stat, 2 * p, lower.tail = FALSE)
  mjb.result = mjb.m(fpca.scores)
  dh1.result = normality.test1(fpca.scores)
  dh2.result = normality.test2(fpca.scores)
  rm(y.fd) # Remove global obejects
  rm(x.fd)
  return(list(p.value.jb, mjb.result$p.value, dh1.result[1], dh2.result[1], p = count))
}

# For regression errors - final function to calculate size (for distr = 'normal') and power of tests
sym.regression.results = function(K = 75, distr = 'normal', counter = 1000, percent = 0.85, param = TRUE)
  # counter - number of replications
  # percent - percent of variation to explain
  # param - T/F; TRUE means 5.1 process; FALSE means 5.2 process
{
  sym.50 = replicate(counter, sym.regression.sk(K = K, distr = distr, percent = percent, param = param))
  a1 = c(sum(sym.50[1, ] < 0.1), sum(sym.50[1, ] < 0.05), sum(sym.50[1, ] < 0.01)) # JB test
  a2 = c(sum(sym.50[2, ] < 0.1), sum(sym.50[2, ] < 0.05), sum(sym.50[2, ] < 0.01)) # MJB_M test
  a3 = c(sum(sym.50[3, ] < 0.1), sum(sym.50[3, ] < 0.05), sum(sym.50[3, ] < 0.01)) # DH1 test
  a4 = c(sum(sym.50[4, ] < 0.1), sum(sym.50[4, ] < 0.05), sum(sym.50[4, ] < 0.01)) # DH2 test
  result.sym.50 = rbind(a1, a2, a3, a4)
  sym.150 = replicate(counter, sym.regression.sk(K = K, distr = distr, N = 150, percent = percent, param = param))
  b1 = c(sum(sym.150[1, ] < 0.1), sum(sym.150[1, ] < 0.05), sum(sym.150[1, ] < 0.01)) # JB test
  b2 = c(sum(sym.150[2,] < 0.1), sum(sym.150[2, ] < 0.05), sum(sym.150[2, ] < 0.01)) # MJB_M test
  b3 = c(sum(sym.150[3,] < 0.1), sum(sym.150[3, ] < 0.05), sum(sym.150[3, ] < 0.01)) # DH1 test
  b4 = c(sum(sym.150[4,] < 0.1), sum(sym.150[4, ] < 0.05), sum(sym.150[4, ] < 0.01)) # DH2 test
  result.sym.150 = rbind(b1, b2, b3, b4)
  sym.450 = replicate(counter, sym.regression.sk(K = K, distr = distr, N = 450, percent = percent, param = param))
  c1 = c(sum(sym.450[1, ] < 0.1), sum(sym.450[1, ] < 0.05), sum(sym.450[1, ] < 0.01)) # JB test
  c2 = c(sum(sym.450[2, ] < 0.1), sum(sym.450[2, ] < 0.05), sum(sym.450[2, ] < 0.01)) # MJB_M test
  c3 = c(sum(sym.450[3, ] < 0.1), sum(sym.450[3, ] < 0.05), sum(sym.450[3, ] < 0.01)) # DH1 test
  c4 = c(sum(sym.450[4, ] < 0.1), sum(sym.450[4, ] < 0.05), sum(sym.450[4, ] < 0.01)) # DH2 test
  result.sym.450 = rbind(c1, c2, c3, c4)
  result.sym.all = cbind(result.sym.50, result.sym.150, result.sym.450)
  row.names(result.sym.all) = c('Jarque-Bera',  'MJB_M', 'Doornik-Hansen', 'Lobato-Velasco')
  standard.error <- sqrt(result.sym.all / counter * (1 - result.sym.all / counter) / counter)
  p <- list(unlist(sym.50[5, ]), unlist(sym.150[5, ]), unlist(sym.450[5, ]))
  return(list(round(result.sym.all / counter, 3), round(standard.error, 3), p))
}

# For real data
real.data.normality.test <- function(data, percent = 0.85)
  # data - real data set (observations in rows)
  # percent - percent of variation to explain
{
  N <- nrow(data)
  K <- ncol(data)
  data.time <- seq(0, 1, length = K)
  data.range <- c(0, 1)
  data.basis <- create.bspline.basis(data.range, nbasis = 10)
  data.fd <- smooth.basisPar(data.time, t(data), data.basis)$fd
  fpca.model <- pca.fd(data.fd, nharm = 10) # FPCA model
  count <- min(which(cumsum(fpca.model$varprop) > percent)) # Dimension of functional data projection
  count <- ifelse(count == 1, 2, count)
  fpca.scores <- fpca.model$scores[, 1:count] # FPCA scores
  data.kurt <- sapply(data.frame(fpca.scores), kurtosis)
  data.skewness <- sapply(data.frame(fpca.scores), skewness)
  stat <- N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p <- length(data.skewness)
  jb.result <- pchisq(stat, 2 * p, lower.tail = FALSE)
  mjb.result <- mjb.m(fpca.scores)
  dh1.result <- normality.test1(fpca.scores)
  dh2.result <- normality.test2(fpca.scores)
  return(list(JB = jb.result, MJB = mjb.result$p.value, DH = dh1.result[1], LV = dh2.result[1], p = count))
}

# For real data (sparse)
real.data.normality.test.sparse <- function(fpca.scores)
  # fpca.scores - fpca scores from fdapace::FPCA
{
  data.kurt <- sapply(data.frame(fpca.scores), kurtosis)
  data.skewness <- sapply(data.frame(fpca.scores), skewness)
  N <- nrow(fpca.scores)
  stat <- N * sum((data.skewness^2 / 6 + (data.kurt - 3)^2 / 24))
  p <- length(data.skewness)
  jb.result <- pchisq(stat, 2 * p, lower.tail = FALSE)
  mjb.result <- mjb.m(fpca.scores)
  dh1.result <- normality.test1(fpca.scores)
  dh2.result <- normality.test2(fpca.scores)
  return(list(JB = jb.result, MJB = mjb.result$p.value, DH = dh1.result[1], LV = dh2.result[1], p = p))
}
