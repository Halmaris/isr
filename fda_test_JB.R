source('fda_JB_functions.R')

gaittime = as.numeric(dimnames(gait)[[1]]) * 20
gaitrange = c(0, 20)
gaitbasis = create.fourier.basis(gaitrange, nbasis = 21)
harmaccelLfd = vec2Lfd(c(0, (2 * pi / 20)^2, 0), rangeval = gaitrange)
gaitfd = smooth.basisPar(gaittime, gait, gaitbasis, Lfdobj = harmaccelLfd, lambda = 1e-2)$fd
hipfd  = gaitfd[, 1]
kneefd = gaitfd[, 2]
model.regression.knee = fRegress(kneefd ~ hipfd) # Funcional regression model
model.predict = predict(model.regression.knee) # Prediction from functional regression model

model.residuals = model.predict$fd - kneefd # Matrix of coefficients for residuals (functional object)

plot(kneefd, main = 'Data')
plot(model.predict, main = 'Prediction')
plot(model.residuals, main = 'Residuals')

fpca = pca.fd(model.residuals, nharm = 4) # Functional PCA, nharm - number of components
fpca$harmonics # Eeight curves
fpca$scores # Scores

# Splines
gaittime = as.numeric(dimnames(gait)[[1]]) * 20
gaitrange = c(0, 20)
gaitbasis = create.bspline.basis(gaitrange)
gaitfd = smooth.basisPar(gaittime, gait, gaitbasis)$fd
hipfd  = gaitfd[, 1]
kneefd = gaitfd[, 2]
model.regression.knee = fRegress(kneefd ~ hipfd) # Funcional regression model
model.predict = predict(model.regression.knee) # Prediction from functional regression model

model.residuals = model.predict$fd - kneefd # Matrix of coefficients for residuals (functional object)

plot(kneefd, main = 'Data')
plot(model.predict, main = 'Prediction')
plot(model.residuals, main = 'Residuals')

fpca = pca.fd(model.residuals, nharm = 4) # Functional PCA, nharm - number of components
fpca$harmonics # Weight curves
fpca$scores # Scores

par(mfrow = c(1, 2))
for (i in 1:4)
{
  acf(fpca$scores[, i])
  acf(fpca$scores[, i]^2)
  model.lm = lm(fpca$scores[, i]~1)
  print(durbinWatsonTest(model.lm))
  if (i < 4)
    readline("Press return to continue")
}

mvn(fpca$scores, mvnTest = 'hz', desc = FALSE) # Henze-Zirkler's Multivariate Normality Test
mvn(fpca$scores, mvnTest = 'mardia', desc = FALSE) # Mardia's Multivariate Normality Test
mvn(fpca$scores, mvnTest = 'royston', desc = FALSE) # Royston's Multivariate Normality Test
mvn(fpca$scores, mvnTest = 'energy', desc = FALSE) # Energy test for MVN
mvn(fpca$scores, mvnTest = 'dh', desc = FALSE) # Doornik-Hansen Omnibus Test for Normality with Independence (Jarque Bera test)
normality.test2(fpca$scores) # Doornik-Hansen Omnibus Test for Normality with allowance for the variable(s) being weakly dependent
mjb.m(fpca$scores) # Test Jarque–Bera (MJB_M, MJB_M*)

###################################
# Simulation studies for Random Walks

# Sizes
sym.size.50 = replicate(1000, sym.size()) #sum(sym.size.50<0.1) - sizes of test
a1 = c(sum(sym.size.50[1,]<0.1), sum(sym.size.50[1,]<0.05), sum(sym.size.50[1,]<0.01)) # sizes of single test
a2 = c(sum(sym.size.50[2,]<0.1), sum(sym.size.50[2,]<0.05), sum(sym.size.50[2,]<0.01)) # sizes of HZ test
a3 = c(sum(sym.size.50[3,]<0.1), sum(sym.size.50[3,]<0.05), sum(sym.size.50[3,]<0.01)) # sizes of Mardia skew test
a4 = c(sum(sym.size.50[4,]<0.1), sum(sym.size.50[4,]<0.05), sum(sym.size.50[4,]<0.01)) # sizes of Mardia kurt test
a5 = c(sum(sym.size.50[5,]<0.1), sum(sym.size.50[5,]<0.05), sum(sym.size.50[5,]<0.01)) # sizes of Royston test
a6 = c(sum(sym.size.50[6,]<0.1), sum(sym.size.50[6,]<0.05), sum(sym.size.50[6,]<0.01)) # sizes of MJB_M test
a7 = c(sum(sym.size.50[7,]<0.1), sum(sym.size.50[7,]<0.05), sum(sym.size.50[7,]<0.01)) # sizes of MJB_M* test
a8 = c(sum(sym.size.50[8,]<0.1), sum(sym.size.50[8,]<0.05), sum(sym.size.50[8,]<0.01)) # sizes of energy test
a9 = c(sum(sym.size.50[9,]<0.1), sum(sym.size.50[9,]<0.05), sum(sym.size.50[9,]<0.01)) # sizes of DH1 test
a10 = c(sum(sym.size.50[10,]<0.1), sum(sym.size.50[10,]<0.05), sum(sym.size.50[10,]<0.01)) # sizes of DH2 test
result.sym.size.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.size.150 = replicate(1000, sym.size(N = 150))
b1 = c(sum(sym.size.150[1,]<0.1), sum(sym.size.150[1,]<0.05), sum(sym.size.150[1,]<0.01)) # sizes of single test
b2 = c(sum(sym.size.150[2,]<0.1), sum(sym.size.150[2,]<0.05), sum(sym.size.150[2,]<0.01)) # sizes of HZ test
b3 = c(sum(sym.size.150[3,]<0.1), sum(sym.size.150[3,]<0.05), sum(sym.size.150[3,]<0.01)) # sizes of Mardia skew test
b4 = c(sum(sym.size.150[4,]<0.1), sum(sym.size.150[4,]<0.05), sum(sym.size.150[4,]<0.01)) # sizes of Mardia kurt test
b5 = c(sum(sym.size.150[5,]<0.1), sum(sym.size.150[5,]<0.05), sum(sym.size.150[5,]<0.01)) # sizes of Royston test
b6 = c(sum(sym.size.150[6,]<0.1), sum(sym.size.150[6,]<0.05), sum(sym.size.150[6,]<0.01)) # sizes of MJB_M test
b7 = c(sum(sym.size.150[7,]<0.1), sum(sym.size.150[7,]<0.05), sum(sym.size.150[7,]<0.01)) # sizes of MJB_M* test
b8 = c(sum(sym.size.150[8,]<0.1), sum(sym.size.150[8,]<0.05), sum(sym.size.150[8,]<0.01)) # sizes of energy test
b9 = c(sum(sym.size.150[9,]<0.1), sum(sym.size.150[9,]<0.05), sum(sym.size.150[9,]<0.01)) # sizes of DH1 test
b10 = c(sum(sym.size.150[10,]<0.1), sum(sym.size.150[10,]<0.05), sum(sym.size.150[10,]<0.01)) # sizes of DH2 test
result.sym.size.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.size.450 = replicate(1000, sym.size(N = 450))
c1 = c(sum(sym.size.450[1,]<0.1), sum(sym.size.450[1,]<0.05), sum(sym.size.450[1,]<0.01)) # sizes of single test
c2 = c(sum(sym.size.450[2,]<0.1), sum(sym.size.450[2,]<0.05), sum(sym.size.450[2,]<0.01)) # sizes of HZ test
c3 = c(sum(sym.size.450[3,]<0.1), sum(sym.size.450[3,]<0.05), sum(sym.size.450[3,]<0.01)) # sizes of Mardia skew test
c4 = c(sum(sym.size.450[4,]<0.1), sum(sym.size.450[4,]<0.05), sum(sym.size.450[4,]<0.01)) # sizes of Mardia kurt test
c5 = c(sum(sym.size.450[5,]<0.1), sum(sym.size.450[5,]<0.05), sum(sym.size.450[5,]<0.01)) # sizes of Royston test
c6 = c(sum(sym.size.450[6,]<0.1), sum(sym.size.450[6,]<0.05), sum(sym.size.450[6,]<0.01)) # sizes of MJB_M test
c7 = c(sum(sym.size.450[7,]<0.1), sum(sym.size.450[7,]<0.05), sum(sym.size.450[7,]<0.01)) # sizes of MJB_M* test
c8 = c(sum(sym.size.450[8,]<0.1), sum(sym.size.450[8,]<0.05), sum(sym.size.450[8,]<0.01)) # sizes of energy test
c9 = c(sum(sym.size.450[9,]<0.1), sum(sym.size.450[9,]<0.05), sum(sym.size.450[9,]<0.01)) # sizes of DH1 test
c10 = c(sum(sym.size.450[10,]<0.1), sum(sym.size.450[10,]<0.05), sum(sym.size.450[10,]<0.01)) # sizes of DH2 test
result.sym.size.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.size = cbind(result.sym.size.50, result.sym.size.150, result.sym.size.450)
row.names(result.sym.size) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                               'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.size, numeric.dollar = F) # Export to latex (in working directory)

# Power against t(2)
sym.power.50 = replicate(1000, sym.power()) #sum(sym.power.50<0.1) - power of test
a1 = c(sum(sym.power.50[1,]<0.1), sum(sym.power.50[1,]<0.05), sum(sym.power.50[1,]<0.01)) #powers of single test
a2 = c(sum(sym.power.50[2,]<0.1), sum(sym.power.50[2,]<0.05), sum(sym.power.50[2,]<0.01)) #powers of HZ test
a3 = c(sum(sym.power.50[3,]<0.1), sum(sym.power.50[3,]<0.05), sum(sym.power.50[3,]<0.01)) #powers of Mardia skew test
a4 = c(sum(sym.power.50[4,]<0.1), sum(sym.power.50[4,]<0.05), sum(sym.power.50[4,]<0.01)) #powers of Mardia kurt test
a5 = c(sum(sym.power.50[5,]<0.1), sum(sym.power.50[5,]<0.05), sum(sym.power.50[5,]<0.01)) #powers of Royston test
a6 = c(sum(sym.power.50[6,]<0.1), sum(sym.power.50[6,]<0.05), sum(sym.power.50[6,]<0.01)) #powers of MJB_M test
a7 = c(sum(sym.power.50[7,]<0.1), sum(sym.power.50[7,]<0.05), sum(sym.power.50[7,]<0.01)) #powers of MJB_M* test
a8 = c(sum(sym.power.50[8,]<0.1), sum(sym.power.50[8,]<0.05), sum(sym.power.50[8,]<0.01)) #powers of energy test
a9 = c(sum(sym.power.50[9,]<0.1), sum(sym.power.50[9,]<0.05), sum(sym.power.50[9,]<0.01)) #powers of DH1 test
a10 = c(sum(sym.power.50[10,]<0.1), sum(sym.power.50[10,]<0.05), sum(sym.power.50[10,]<0.01)) #powers of DH2 test
result.sym.power.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.power.150 = replicate(1000, sym.power(N = 150))
b1 = c(sum(sym.power.150[1,]<0.1), sum(sym.power.150[1,]<0.05), sum(sym.power.150[1,]<0.01)) #powers of single test
b2 = c(sum(sym.power.150[2,]<0.1), sum(sym.power.150[2,]<0.05), sum(sym.power.150[2,]<0.01)) #powers of HZ test
b3 = c(sum(sym.power.150[3,]<0.1), sum(sym.power.150[3,]<0.05), sum(sym.power.150[3,]<0.01)) #powers of Mardia skew test
b4 = c(sum(sym.power.150[4,]<0.1), sum(sym.power.150[4,]<0.05), sum(sym.power.150[4,]<0.01)) #powers of Mardia kurt test
b5 = c(sum(sym.power.150[5,]<0.1), sum(sym.power.150[5,]<0.05), sum(sym.power.150[5,]<0.01)) #powers of Royston test
b6 = c(sum(sym.power.150[6,]<0.1), sum(sym.power.150[6,]<0.05), sum(sym.power.150[6,]<0.01)) #powers of MJB_M test
b7 = c(sum(sym.power.150[7,]<0.1), sum(sym.power.150[7,]<0.05), sum(sym.power.150[7,]<0.01)) #powers of MJB_M* test
b8 = c(sum(sym.power.150[8,]<0.1), sum(sym.power.150[8,]<0.05), sum(sym.power.150[8,]<0.01)) #powers of energy test
b9 = c(sum(sym.power.150[9,]<0.1), sum(sym.power.150[9,]<0.05), sum(sym.power.150[9,]<0.01)) #powers of DH1 test
b10 = c(sum(sym.power.150[10,]<0.1), sum(sym.power.150[10,]<0.05), sum(sym.power.150[10,]<0.01)) #powers of DH2 test
result.sym.power.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.power.450 = replicate(1000, sym.power(N = 450))
c1 = c(sum(sym.power.450[1,]<0.1), sum(sym.power.450[1,]<0.05), sum(sym.power.450[1,]<0.01)) #powers of single test
c2 = c(sum(sym.power.450[2,]<0.1), sum(sym.power.450[2,]<0.05), sum(sym.power.450[2,]<0.01)) #powers of HZ test
c3 = c(sum(sym.power.450[3,]<0.1), sum(sym.power.450[3,]<0.05), sum(sym.power.450[3,]<0.01)) #powers of Mardia skew test
c4 = c(sum(sym.power.450[4,]<0.1), sum(sym.power.450[4,]<0.05), sum(sym.power.450[4,]<0.01)) #powers of Mardia kurt test
c5 = c(sum(sym.power.450[5,]<0.1), sum(sym.power.450[5,]<0.05), sum(sym.power.450[5,]<0.01)) #powers of Royston test
c6 = c(sum(sym.power.450[6,]<0.1), sum(sym.power.450[6,]<0.05), sum(sym.power.450[6,]<0.01)) #powers of MJB_M test
c7 = c(sum(sym.power.450[7,]<0.1), sum(sym.power.450[7,]<0.05), sum(sym.power.450[7,]<0.01)) #powers of MJB_M* test
c8 = c(sum(sym.power.450[8,]<0.1), sum(sym.power.450[8,]<0.05), sum(sym.power.450[8,]<0.01)) #powers of energy test
c9 = c(sum(sym.power.450[9,]<0.1), sum(sym.power.450[9,]<0.05), sum(sym.power.450[9,]<0.01)) #powers of DH1 test
c10 = c(sum(sym.power.450[10,]<0.1), sum(sym.power.450[10,]<0.05), sum(sym.power.450[10,]<0.01)) #powers of DH2 test
result.sym.power.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.power = cbind(result.sym.power.50, result.sym.power.150, result.sym.power.450)
row.names(result.sym.power) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.power, numeric.dollar = F) #export to latex (in working directory)

# Power against t(4)
sym.power.50 = replicate(1000, sym.power(df = 4)) #sum(sym.power.50<0.1) - power of test
a1 = c(sum(sym.power.50[1,]<0.1), sum(sym.power.50[1,]<0.05), sum(sym.power.50[1,]<0.01)) #powers of single test
a2 = c(sum(sym.power.50[2,]<0.1), sum(sym.power.50[2,]<0.05), sum(sym.power.50[2,]<0.01)) #powers of HZ test
a3 = c(sum(sym.power.50[3,]<0.1), sum(sym.power.50[3,]<0.05), sum(sym.power.50[3,]<0.01)) #powers of Mardia skew test
a4 = c(sum(sym.power.50[4,]<0.1), sum(sym.power.50[4,]<0.05), sum(sym.power.50[4,]<0.01)) #powers of Mardia kurt test
a5 = c(sum(sym.power.50[5,]<0.1), sum(sym.power.50[5,]<0.05), sum(sym.power.50[5,]<0.01)) #powers of Royston test
a6 = c(sum(sym.power.50[6,]<0.1), sum(sym.power.50[6,]<0.05), sum(sym.power.50[6,]<0.01)) #powers of MJB_M test
a7 = c(sum(sym.power.50[7,]<0.1), sum(sym.power.50[7,]<0.05), sum(sym.power.50[7,]<0.01)) #powers of MJB_M* test
a8 = c(sum(sym.power.50[8,]<0.1), sum(sym.power.50[8,]<0.05), sum(sym.power.50[8,]<0.01)) #powers of energy test
a9 = c(sum(sym.power.50[9,]<0.1), sum(sym.power.50[9,]<0.05), sum(sym.power.50[9,]<0.01)) #powers of DH1 test
a10 = c(sum(sym.power.50[10,]<0.1), sum(sym.power.50[10,]<0.05), sum(sym.power.50[10,]<0.01)) #powers of DH2 test
result.sym.power.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.power.150 = replicate(1000, sym.power(N = 150, df = 4))
b1 = c(sum(sym.power.150[1,]<0.1), sum(sym.power.150[1,]<0.05), sum(sym.power.150[1,]<0.01)) #powers of single test
b2 = c(sum(sym.power.150[2,]<0.1), sum(sym.power.150[2,]<0.05), sum(sym.power.150[2,]<0.01)) #powers of HZ test
b3 = c(sum(sym.power.150[3,]<0.1), sum(sym.power.150[3,]<0.05), sum(sym.power.150[3,]<0.01)) #powers of Mardia skew test
b4 = c(sum(sym.power.150[4,]<0.1), sum(sym.power.150[4,]<0.05), sum(sym.power.150[4,]<0.01)) #powers of Mardia kurt test
b5 = c(sum(sym.power.150[5,]<0.1), sum(sym.power.150[5,]<0.05), sum(sym.power.150[5,]<0.01)) #powers of Royston test
b6 = c(sum(sym.power.150[6,]<0.1), sum(sym.power.150[6,]<0.05), sum(sym.power.150[6,]<0.01)) #powers of MJB_M test
b7 = c(sum(sym.power.150[7,]<0.1), sum(sym.power.150[7,]<0.05), sum(sym.power.150[7,]<0.01)) #powers of MJB_M* test
b8 = c(sum(sym.power.150[8,]<0.1), sum(sym.power.150[8,]<0.05), sum(sym.power.150[8,]<0.01)) #powers of energy test
b9 = c(sum(sym.power.150[9,]<0.1), sum(sym.power.150[9,]<0.05), sum(sym.power.150[9,]<0.01)) #powers of DH1 test
b10 = c(sum(sym.power.150[10,]<0.1), sum(sym.power.150[10,]<0.05), sum(sym.power.150[10,]<0.01)) #powers of DH2 test
result.sym.power.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.power.450 = replicate(1000, sym.power(N = 450, df = 4))
c1 = c(sum(sym.power.450[1,]<0.1), sum(sym.power.450[1,]<0.05), sum(sym.power.450[1,]<0.01)) #powers of single test
c2 = c(sum(sym.power.450[2,]<0.1), sum(sym.power.450[2,]<0.05), sum(sym.power.450[2,]<0.01)) #powers of HZ test
c3 = c(sum(sym.power.450[3,]<0.1), sum(sym.power.450[3,]<0.05), sum(sym.power.450[3,]<0.01)) #powers of Mardia skew test
c4 = c(sum(sym.power.450[4,]<0.1), sum(sym.power.450[4,]<0.05), sum(sym.power.450[4,]<0.01)) #powers of Mardia kurt test
c5 = c(sum(sym.power.450[5,]<0.1), sum(sym.power.450[5,]<0.05), sum(sym.power.450[5,]<0.01)) #powers of Royston test
c6 = c(sum(sym.power.450[6,]<0.1), sum(sym.power.450[6,]<0.05), sum(sym.power.450[6,]<0.01)) #powers of MJB_M test
c7 = c(sum(sym.power.450[7,]<0.1), sum(sym.power.450[7,]<0.05), sum(sym.power.450[7,]<0.01)) #powers of MJB_M* test
c8 = c(sum(sym.power.450[8,]<0.1), sum(sym.power.450[8,]<0.05), sum(sym.power.450[8,]<0.01)) #powers of energy test
c9 = c(sum(sym.power.450[9,]<0.1), sum(sym.power.450[9,]<0.05), sum(sym.power.450[9,]<0.01)) #powers of DH1 test
c10 = c(sum(sym.power.450[10,]<0.1), sum(sym.power.450[10,]<0.05), sum(sym.power.450[10,]<0.01)) #powers of DH2 test
result.sym.power.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.power = cbind(result.sym.power.50, result.sym.power.150, result.sym.power.450)
row.names(result.sym.power) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.power, numeric.dollar = F) #export to latex (in working directory)

# Power against t(8)
sym.power.50 = replicate(1000, sym.power(df = 8)) #sum(sym.power.50<0.1) - power of test
a1 = c(sum(sym.power.50[1,]<0.1), sum(sym.power.50[1,]<0.05), sum(sym.power.50[1,]<0.01)) #powers of single test
a2 = c(sum(sym.power.50[2,]<0.1), sum(sym.power.50[2,]<0.05), sum(sym.power.50[2,]<0.01)) #powers of HZ test
a3 = c(sum(sym.power.50[3,]<0.1), sum(sym.power.50[3,]<0.05), sum(sym.power.50[3,]<0.01)) #powers of Mardia skew test
a4 = c(sum(sym.power.50[4,]<0.1), sum(sym.power.50[4,]<0.05), sum(sym.power.50[4,]<0.01)) #powers of Mardia kurt test
a5 = c(sum(sym.power.50[5,]<0.1), sum(sym.power.50[5,]<0.05), sum(sym.power.50[5,]<0.01)) #powers of Royston test
a6 = c(sum(sym.power.50[6,]<0.1), sum(sym.power.50[6,]<0.05), sum(sym.power.50[6,]<0.01)) #powers of MJB_M test
a7 = c(sum(sym.power.50[7,]<0.1), sum(sym.power.50[7,]<0.05), sum(sym.power.50[7,]<0.01)) #powers of MJB_M* test
a8 = c(sum(sym.power.50[8,]<0.1), sum(sym.power.50[8,]<0.05), sum(sym.power.50[8,]<0.01)) #powers of energy test
a9 = c(sum(sym.power.50[9,]<0.1), sum(sym.power.50[9,]<0.05), sum(sym.power.50[9,]<0.01)) #powers of DH1 test
a10 = c(sum(sym.power.50[10,]<0.1), sum(sym.power.50[10,]<0.05), sum(sym.power.50[10,]<0.01)) #powers of DH2 test
result.sym.power.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.power.150 = replicate(1000, sym.power(N = 150, df = 8))
b1 = c(sum(sym.power.150[1,]<0.1), sum(sym.power.150[1,]<0.05), sum(sym.power.150[1,]<0.01)) #powers of single test
b2 = c(sum(sym.power.150[2,]<0.1), sum(sym.power.150[2,]<0.05), sum(sym.power.150[2,]<0.01)) #powers of HZ test
b3 = c(sum(sym.power.150[3,]<0.1), sum(sym.power.150[3,]<0.05), sum(sym.power.150[3,]<0.01)) #powers of Mardia skew test
b4 = c(sum(sym.power.150[4,]<0.1), sum(sym.power.150[4,]<0.05), sum(sym.power.150[4,]<0.01)) #powers of Mardia kurt test
b5 = c(sum(sym.power.150[5,]<0.1), sum(sym.power.150[5,]<0.05), sum(sym.power.150[5,]<0.01)) #powers of Royston test
b6 = c(sum(sym.power.150[6,]<0.1), sum(sym.power.150[6,]<0.05), sum(sym.power.150[6,]<0.01)) #powers of MJB_M test
b7 = c(sum(sym.power.150[7,]<0.1), sum(sym.power.150[7,]<0.05), sum(sym.power.150[7,]<0.01)) #powers of MJB_M* test
b8 = c(sum(sym.power.150[8,]<0.1), sum(sym.power.150[8,]<0.05), sum(sym.power.150[8,]<0.01)) #powers of energy test
b9 = c(sum(sym.power.150[9,]<0.1), sum(sym.power.150[9,]<0.05), sum(sym.power.150[9,]<0.01)) #powers of DH1 test
b10 = c(sum(sym.power.150[10,]<0.1), sum(sym.power.150[10,]<0.05), sum(sym.power.150[10,]<0.01)) #powers of DH2 test
result.sym.power.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.power.450 = replicate(1000, sym.power(N = 450, df = 8))
c1 = c(sum(sym.power.450[1,]<0.1), sum(sym.power.450[1,]<0.05), sum(sym.power.450[1,]<0.01)) #powers of single test
c2 = c(sum(sym.power.450[2,]<0.1), sum(sym.power.450[2,]<0.05), sum(sym.power.450[2,]<0.01)) #powers of HZ test
c3 = c(sum(sym.power.450[3,]<0.1), sum(sym.power.450[3,]<0.05), sum(sym.power.450[3,]<0.01)) #powers of Mardia skew test
c4 = c(sum(sym.power.450[4,]<0.1), sum(sym.power.450[4,]<0.05), sum(sym.power.450[4,]<0.01)) #powers of Mardia kurt test
c5 = c(sum(sym.power.450[5,]<0.1), sum(sym.power.450[5,]<0.05), sum(sym.power.450[5,]<0.01)) #powers of Royston test
c6 = c(sum(sym.power.450[6,]<0.1), sum(sym.power.450[6,]<0.05), sum(sym.power.450[6,]<0.01)) #powers of MJB_M test
c7 = c(sum(sym.power.450[7,]<0.1), sum(sym.power.450[7,]<0.05), sum(sym.power.450[7,]<0.01)) #powers of MJB_M* test
c8 = c(sum(sym.power.450[8,]<0.1), sum(sym.power.450[8,]<0.05), sum(sym.power.450[8,]<0.01)) #powers of energy test
c9 = c(sum(sym.power.450[9,]<0.1), sum(sym.power.450[9,]<0.05), sum(sym.power.450[9,]<0.01)) #powers of DH1 test
c10 = c(sum(sym.power.450[10,]<0.1), sum(sym.power.450[10,]<0.05), sum(sym.power.450[10,]<0.01)) #powers of DH2 test
result.sym.power.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.power = cbind(result.sym.power.50, result.sym.power.150, result.sym.power.450)
row.names(result.sym.power) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.power, numeric.dollar = F) # Export to latex (in working directory)

###################################
# Symulations studies for regression

# Sizes (regression)
sym.size.reg.50 = replicate(1000, sym.size.reg())
a1 = c(sum(sym.size.reg.50[1,]<0.1), sum(sym.size.reg.50[1,]<0.05), sum(sym.size.reg.50[1,]<0.01)) #sizes of single test
a2 = c(sum(sym.size.reg.50[2,]<0.1), sum(sym.size.reg.50[2,]<0.05), sum(sym.size.reg.50[2,]<0.01)) #sizes of HZ test
a3 = c(sum(sym.size.reg.50[3,]<0.1), sum(sym.size.reg.50[3,]<0.05), sum(sym.size.reg.50[3,]<0.01)) #sizes of Mardia skew test
a4 = c(sum(sym.size.reg.50[4,]<0.1), sum(sym.size.reg.50[4,]<0.05), sum(sym.size.reg.50[4,]<0.01)) #sizes of Mardia kurt test
a5 = c(sum(sym.size.reg.50[5,]<0.1), sum(sym.size.reg.50[5,]<0.05), sum(sym.size.reg.50[5,]<0.01)) #sizes of Royston test
a6 = c(sum(sym.size.reg.50[6,]<0.1), sum(sym.size.reg.50[6,]<0.05), sum(sym.size.reg.50[6,]<0.01)) #sizes of MJB_M test
a7 = c(sum(sym.size.reg.50[7,]<0.1), sum(sym.size.reg.50[7,]<0.05), sum(sym.size.reg.50[7,]<0.01)) #sizes of MJB_M* test
a8 = c(sum(sym.size.reg.50[8,]<0.1), sum(sym.size.reg.50[8,]<0.05), sum(sym.size.reg.50[8,]<0.01)) #sizes of energy test
a9 = c(sum(sym.size.reg.50[9,]<0.1), sum(sym.size.reg.50[9,]<0.05), sum(sym.size.reg.50[9,]<0.01)) #sizes of DH1 test
a10 = c(sum(sym.size.reg.50[10,]<0.1), sum(sym.size.reg.50[10,]<0.05), sum(sym.size.reg.50[10,]<0.01)) #sizes of DH2 test
result.sym.size.reg.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.size.reg.150 = replicate(1000, sym.size.reg(N = 150))
b1 = c(sum(sym.size.reg.150[1,]<0.1), sum(sym.size.reg.150[1,]<0.05), sum(sym.size.reg.150[1,]<0.01)) #sizes of single test
b2 = c(sum(sym.size.reg.150[2,]<0.1), sum(sym.size.reg.150[2,]<0.05), sum(sym.size.reg.150[2,]<0.01)) #sizes of HZ test
b3 = c(sum(sym.size.reg.150[3,]<0.1), sum(sym.size.reg.150[3,]<0.05), sum(sym.size.reg.150[3,]<0.01)) #sizes of Mardia skew test
b4 = c(sum(sym.size.reg.150[4,]<0.1), sum(sym.size.reg.150[4,]<0.05), sum(sym.size.reg.150[4,]<0.01)) #sizes of Mardia kurt test
b5 = c(sum(sym.size.reg.150[5,]<0.1), sum(sym.size.reg.150[5,]<0.05), sum(sym.size.reg.150[5,]<0.01)) #sizes of Royston test
b6 = c(sum(sym.size.reg.150[6,]<0.1), sum(sym.size.reg.150[6,]<0.05), sum(sym.size.reg.150[6,]<0.01)) #sizes of MJB_M test
b7 = c(sum(sym.size.reg.150[7,]<0.1), sum(sym.size.reg.150[7,]<0.05), sum(sym.size.reg.150[7,]<0.01)) #sizes of MJB_M* test
b8 = c(sum(sym.size.reg.150[8,]<0.1), sum(sym.size.reg.150[8,]<0.05), sum(sym.size.reg.150[8,]<0.01)) #sizes of energy test
b9 = c(sum(sym.size.reg.150[9,]<0.1), sum(sym.size.reg.150[9,]<0.05), sum(sym.size.reg.150[9,]<0.01)) #sizes of DH1 test
b10 = c(sum(sym.size.reg.150[10,]<0.1), sum(sym.size.reg.150[10,]<0.05), sum(sym.size.reg.150[10,]<0.01)) #sizes of DH2 test
result.sym.size.reg.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.size.reg.450 = replicate(1000, sym.size.reg(N = 450))
c1 = c(sum(sym.size.reg.450[1,]<0.1), sum(sym.size.reg.450[1,]<0.05), sum(sym.size.reg.450[1,]<0.01)) #sizes of single test
c2 = c(sum(sym.size.reg.450[2,]<0.1), sum(sym.size.reg.450[2,]<0.05), sum(sym.size.reg.450[2,]<0.01)) #sizes of HZ test
c3 = c(sum(sym.size.reg.450[3,]<0.1), sum(sym.size.reg.450[3,]<0.05), sum(sym.size.reg.450[3,]<0.01)) #sizes of Mardia skew test
c4 = c(sum(sym.size.reg.450[4,]<0.1), sum(sym.size.reg.450[4,]<0.05), sum(sym.size.reg.450[4,]<0.01)) #sizes of Mardia kurt test
c5 = c(sum(sym.size.reg.450[5,]<0.1), sum(sym.size.reg.450[5,]<0.05), sum(sym.size.reg.450[5,]<0.01)) #sizes of Royston test
c6 = c(sum(sym.size.reg.450[6,]<0.1), sum(sym.size.reg.450[6,]<0.05), sum(sym.size.reg.450[6,]<0.01)) #sizes of MJB_M test
c7 = c(sum(sym.size.reg.450[7,]<0.1), sum(sym.size.reg.450[7,]<0.05), sum(sym.size.reg.450[7,]<0.01)) #sizes of MJB_M* test
c8 = c(sum(sym.size.reg.450[8,]<0.1), sum(sym.size.reg.450[8,]<0.05), sum(sym.size.reg.450[8,]<0.01)) #sizes of energy test
c9 = c(sum(sym.size.reg.450[9,]<0.1), sum(sym.size.reg.450[9,]<0.05), sum(sym.size.reg.450[9,]<0.01)) #sizes of DH1 test
c10 = c(sum(sym.size.reg.450[10,]<0.1), sum(sym.size.reg.450[10,]<0.05), sum(sym.size.reg.450[10,]<0.01)) #sizes of DH2 test
result.sym.size.reg.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.size.reg = cbind(result.sym.size.reg.50, result.sym.size.reg.150, result.sym.size.reg.450)
row.names(result.sym.size.reg) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                   'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.size.reg, numeric.dollar = F) #export to latex (in working directory)

# Powers (regression) df = 2
sym.size.reg.50 = replicate(1000, sym.size.reg(size = FALSE))
a1 = c(sum(sym.size.reg.50[1,]<0.1), sum(sym.size.reg.50[1,]<0.05), sum(sym.size.reg.50[1,]<0.01)) #powers of single test
a2 = c(sum(sym.size.reg.50[2,]<0.1), sum(sym.size.reg.50[2,]<0.05), sum(sym.size.reg.50[2,]<0.01)) #powers of HZ test
a3 = c(sum(sym.size.reg.50[3,]<0.1), sum(sym.size.reg.50[3,]<0.05), sum(sym.size.reg.50[3,]<0.01)) #powers of Mardia skew test
a4 = c(sum(sym.size.reg.50[4,]<0.1), sum(sym.size.reg.50[4,]<0.05), sum(sym.size.reg.50[4,]<0.01)) #powers of Mardia kurt test
a5 = c(sum(sym.size.reg.50[5,]<0.1), sum(sym.size.reg.50[5,]<0.05), sum(sym.size.reg.50[5,]<0.01)) #powers of Royston test
a6 = c(sum(sym.size.reg.50[6,]<0.1), sum(sym.size.reg.50[6,]<0.05), sum(sym.size.reg.50[6,]<0.01)) #powers of MJB_M test
a7 = c(sum(sym.size.reg.50[7,]<0.1), sum(sym.size.reg.50[7,]<0.05), sum(sym.size.reg.50[7,]<0.01)) #powers of MJB_M* test
a8 = c(sum(sym.size.reg.50[8,]<0.1), sum(sym.size.reg.50[8,]<0.05), sum(sym.size.reg.50[8,]<0.01)) #powers of energy test
a9 = c(sum(sym.size.reg.50[9,]<0.1), sum(sym.size.reg.50[9,]<0.05), sum(sym.size.reg.50[9,]<0.01)) #powers of DH1 test
a10 = c(sum(sym.size.reg.50[10,]<0.1), sum(sym.size.reg.50[10,]<0.05), sum(sym.size.reg.50[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.size.reg.150 = replicate(1000, sym.size.reg(N = 150, size = FALSE))
b1 = c(sum(sym.size.reg.150[1,]<0.1), sum(sym.size.reg.150[1,]<0.05), sum(sym.size.reg.150[1,]<0.01)) #powers of single test
b2 = c(sum(sym.size.reg.150[2,]<0.1), sum(sym.size.reg.150[2,]<0.05), sum(sym.size.reg.150[2,]<0.01)) #powers of HZ test
b3 = c(sum(sym.size.reg.150[3,]<0.1), sum(sym.size.reg.150[3,]<0.05), sum(sym.size.reg.150[3,]<0.01)) #powers of Mardia skew test
b4 = c(sum(sym.size.reg.150[4,]<0.1), sum(sym.size.reg.150[4,]<0.05), sum(sym.size.reg.150[4,]<0.01)) #powers of Mardia kurt test
b5 = c(sum(sym.size.reg.150[5,]<0.1), sum(sym.size.reg.150[5,]<0.05), sum(sym.size.reg.150[5,]<0.01)) #powers of Royston test
b6 = c(sum(sym.size.reg.150[6,]<0.1), sum(sym.size.reg.150[6,]<0.05), sum(sym.size.reg.150[6,]<0.01)) #powers of MJB_M test
b7 = c(sum(sym.size.reg.150[7,]<0.1), sum(sym.size.reg.150[7,]<0.05), sum(sym.size.reg.150[7,]<0.01)) #powers of MJB_M* test
b8 = c(sum(sym.size.reg.150[8,]<0.1), sum(sym.size.reg.150[8,]<0.05), sum(sym.size.reg.150[8,]<0.01)) #powers of energy test
b9 = c(sum(sym.size.reg.150[9,]<0.1), sum(sym.size.reg.150[9,]<0.05), sum(sym.size.reg.150[9,]<0.01)) #powers of DH1 test
b10 = c(sum(sym.size.reg.150[10,]<0.1), sum(sym.size.reg.150[10,]<0.05), sum(sym.size.reg.150[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.size.reg.450 = replicate(1000, sym.size.reg(N = 450, size = FALSE))
c1 = c(sum(sym.size.reg.450[1,]<0.1), sum(sym.size.reg.450[1,]<0.05), sum(sym.size.reg.450[1,]<0.01)) #powers of single test
c2 = c(sum(sym.size.reg.450[2,]<0.1), sum(sym.size.reg.450[2,]<0.05), sum(sym.size.reg.450[2,]<0.01)) #powers of HZ test
c3 = c(sum(sym.size.reg.450[3,]<0.1), sum(sym.size.reg.450[3,]<0.05), sum(sym.size.reg.450[3,]<0.01)) #powers of Mardia skew test
c4 = c(sum(sym.size.reg.450[4,]<0.1), sum(sym.size.reg.450[4,]<0.05), sum(sym.size.reg.450[4,]<0.01)) #powers of Mardia kurt test
c5 = c(sum(sym.size.reg.450[5,]<0.1), sum(sym.size.reg.450[5,]<0.05), sum(sym.size.reg.450[5,]<0.01)) #powers of Royston test
c6 = c(sum(sym.size.reg.450[6,]<0.1), sum(sym.size.reg.450[6,]<0.05), sum(sym.size.reg.450[6,]<0.01)) #powers of MJB_M test
c7 = c(sum(sym.size.reg.450[7,]<0.1), sum(sym.size.reg.450[7,]<0.05), sum(sym.size.reg.450[7,]<0.01)) #powers of MJB_M* test
c8 = c(sum(sym.size.reg.450[8,]<0.1), sum(sym.size.reg.450[8,]<0.05), sum(sym.size.reg.450[8,]<0.01)) #powers of energy test
c9 = c(sum(sym.size.reg.450[9,]<0.1), sum(sym.size.reg.450[9,]<0.05), sum(sym.size.reg.450[9,]<0.01)) #powers of DH1 test
c10 = c(sum(sym.size.reg.450[10,]<0.1), sum(sym.size.reg.450[10,]<0.05), sum(sym.size.reg.450[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.size.reg = cbind(result.sym.size.reg.50, result.sym.size.reg.150, result.sym.size.reg.450)
row.names(result.sym.size.reg) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.size.reg, numeric.dollar = FALSE) #export to latex (in working directory)

# df = 4
sym.size.reg.50 = replicate(1000, sym.size.reg(df = 4, size = FALSE))
a1 = c(sum(sym.size.reg.50[1,]<0.1), sum(sym.size.reg.50[1,]<0.05), sum(sym.size.reg.50[1,]<0.01)) #powers of single test
a2 = c(sum(sym.size.reg.50[2,]<0.1), sum(sym.size.reg.50[2,]<0.05), sum(sym.size.reg.50[2,]<0.01)) #powers of HZ test
a3 = c(sum(sym.size.reg.50[3,]<0.1), sum(sym.size.reg.50[3,]<0.05), sum(sym.size.reg.50[3,]<0.01)) #powers of Mardia skew test
a4 = c(sum(sym.size.reg.50[4,]<0.1), sum(sym.size.reg.50[4,]<0.05), sum(sym.size.reg.50[4,]<0.01)) #powers of Mardia kurt test
a5 = c(sum(sym.size.reg.50[5,]<0.1), sum(sym.size.reg.50[5,]<0.05), sum(sym.size.reg.50[5,]<0.01)) #powers of Royston test
a6 = c(sum(sym.size.reg.50[6,]<0.1), sum(sym.size.reg.50[6,]<0.05), sum(sym.size.reg.50[6,]<0.01)) #powers of MJB_M test
a7 = c(sum(sym.size.reg.50[7,]<0.1), sum(sym.size.reg.50[7,]<0.05), sum(sym.size.reg.50[7,]<0.01)) #powers of MJB_M* test
a8 = c(sum(sym.size.reg.50[8,]<0.1), sum(sym.size.reg.50[8,]<0.05), sum(sym.size.reg.50[8,]<0.01)) #powers of energy test
a9 = c(sum(sym.size.reg.50[9,]<0.1), sum(sym.size.reg.50[9,]<0.05), sum(sym.size.reg.50[9,]<0.01)) #powers of DH1 test
a10 = c(sum(sym.size.reg.50[10,]<0.1), sum(sym.size.reg.50[10,]<0.05), sum(sym.size.reg.50[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.size.reg.150 = replicate(1000, sym.size.reg(N = 150, df = 4, size = FALSE))
b1 = c(sum(sym.size.reg.150[1,]<0.1), sum(sym.size.reg.150[1,]<0.05), sum(sym.size.reg.150[1,]<0.01)) #powers of single test
b2 = c(sum(sym.size.reg.150[2,]<0.1), sum(sym.size.reg.150[2,]<0.05), sum(sym.size.reg.150[2,]<0.01)) #powers of HZ test
b3 = c(sum(sym.size.reg.150[3,]<0.1), sum(sym.size.reg.150[3,]<0.05), sum(sym.size.reg.150[3,]<0.01)) #powers of Mardia skew test
b4 = c(sum(sym.size.reg.150[4,]<0.1), sum(sym.size.reg.150[4,]<0.05), sum(sym.size.reg.150[4,]<0.01)) #powers of Mardia kurt test
b5 = c(sum(sym.size.reg.150[5,]<0.1), sum(sym.size.reg.150[5,]<0.05), sum(sym.size.reg.150[5,]<0.01)) #powers of Royston test
b6 = c(sum(sym.size.reg.150[6,]<0.1), sum(sym.size.reg.150[6,]<0.05), sum(sym.size.reg.150[6,]<0.01)) #powers of MJB_M test
b7 = c(sum(sym.size.reg.150[7,]<0.1), sum(sym.size.reg.150[7,]<0.05), sum(sym.size.reg.150[7,]<0.01)) #powers of MJB_M* test
b8 = c(sum(sym.size.reg.150[8,]<0.1), sum(sym.size.reg.150[8,]<0.05), sum(sym.size.reg.150[8,]<0.01)) #powers of energy test
b9 = c(sum(sym.size.reg.150[9,]<0.1), sum(sym.size.reg.150[9,]<0.05), sum(sym.size.reg.150[9,]<0.01)) #powers of DH1 test
b10 = c(sum(sym.size.reg.150[10,]<0.1), sum(sym.size.reg.150[10,]<0.05), sum(sym.size.reg.150[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.size.reg.450 = replicate(1000, sym.size.reg(N = 450, df = 4, size = FALSE))
c1 = c(sum(sym.size.reg.450[1,]<0.1), sum(sym.size.reg.450[1,]<0.05), sum(sym.size.reg.450[1,]<0.01)) #powers of single test
c2 = c(sum(sym.size.reg.450[2,]<0.1), sum(sym.size.reg.450[2,]<0.05), sum(sym.size.reg.450[2,]<0.01)) #powers of HZ test
c3 = c(sum(sym.size.reg.450[3,]<0.1), sum(sym.size.reg.450[3,]<0.05), sum(sym.size.reg.450[3,]<0.01)) #powers of Mardia skew test
c4 = c(sum(sym.size.reg.450[4,]<0.1), sum(sym.size.reg.450[4,]<0.05), sum(sym.size.reg.450[4,]<0.01)) #powers of Mardia kurt test
c5 = c(sum(sym.size.reg.450[5,]<0.1), sum(sym.size.reg.450[5,]<0.05), sum(sym.size.reg.450[5,]<0.01)) #powers of Royston test
c6 = c(sum(sym.size.reg.450[6,]<0.1), sum(sym.size.reg.450[6,]<0.05), sum(sym.size.reg.450[6,]<0.01)) #powers of MJB_M test
c7 = c(sum(sym.size.reg.450[7,]<0.1), sum(sym.size.reg.450[7,]<0.05), sum(sym.size.reg.450[7,]<0.01)) #powers of MJB_M* test
c8 = c(sum(sym.size.reg.450[8,]<0.1), sum(sym.size.reg.450[8,]<0.05), sum(sym.size.reg.450[8,]<0.01)) #powers of energy test
c9 = c(sum(sym.size.reg.450[9,]<0.1), sum(sym.size.reg.450[9,]<0.05), sum(sym.size.reg.450[9,]<0.01)) #powers of DH1 test
c10 = c(sum(sym.size.reg.450[10,]<0.1), sum(sym.size.reg.450[10,]<0.05), sum(sym.size.reg.450[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.size.reg = cbind(result.sym.size.reg.50, result.sym.size.reg.150, result.sym.size.reg.450)
row.names(result.sym.size.reg) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                   'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.size.reg, numeric.dollar = FALSE) #export to latex (in working directory)

# df = 8
sym.size.reg.50 = replicate(1000, sym.size.reg(df = 8, size = FALSE))
a1 = c(sum(sym.size.reg.50[1,]<0.1), sum(sym.size.reg.50[1,]<0.05), sum(sym.size.reg.50[1,]<0.01)) #powers of single test
a2 = c(sum(sym.size.reg.50[2,]<0.1), sum(sym.size.reg.50[2,]<0.05), sum(sym.size.reg.50[2,]<0.01)) #powers of HZ test
a3 = c(sum(sym.size.reg.50[3,]<0.1), sum(sym.size.reg.50[3,]<0.05), sum(sym.size.reg.50[3,]<0.01)) #powers of Mardia skew test
a4 = c(sum(sym.size.reg.50[4,]<0.1), sum(sym.size.reg.50[4,]<0.05), sum(sym.size.reg.50[4,]<0.01)) #powers of Mardia kurt test
a5 = c(sum(sym.size.reg.50[5,]<0.1), sum(sym.size.reg.50[5,]<0.05), sum(sym.size.reg.50[5,]<0.01)) #powers of Royston test
a6 = c(sum(sym.size.reg.50[6,]<0.1), sum(sym.size.reg.50[6,]<0.05), sum(sym.size.reg.50[6,]<0.01)) #powers of MJB_M test
a7 = c(sum(sym.size.reg.50[7,]<0.1), sum(sym.size.reg.50[7,]<0.05), sum(sym.size.reg.50[7,]<0.01)) #powers of MJB_M* test
a8 = c(sum(sym.size.reg.50[8,]<0.1), sum(sym.size.reg.50[8,]<0.05), sum(sym.size.reg.50[8,]<0.01)) #powers of energy test
a9 = c(sum(sym.size.reg.50[9,]<0.1), sum(sym.size.reg.50[9,]<0.05), sum(sym.size.reg.50[9,]<0.01)) #powers of DH1 test
a10 = c(sum(sym.size.reg.50[10,]<0.1), sum(sym.size.reg.50[10,]<0.05), sum(sym.size.reg.50[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.50 = rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
sym.size.reg.150 = replicate(1000, sym.size.reg(N = 150, df = 8, size = FALSE))
b1 = c(sum(sym.size.reg.150[1,]<0.1), sum(sym.size.reg.150[1,]<0.05), sum(sym.size.reg.150[1,]<0.01)) #powers of single test
b2 = c(sum(sym.size.reg.150[2,]<0.1), sum(sym.size.reg.150[2,]<0.05), sum(sym.size.reg.150[2,]<0.01)) #powers of HZ test
b3 = c(sum(sym.size.reg.150[3,]<0.1), sum(sym.size.reg.150[3,]<0.05), sum(sym.size.reg.150[3,]<0.01)) #powers of Mardia skew test
b4 = c(sum(sym.size.reg.150[4,]<0.1), sum(sym.size.reg.150[4,]<0.05), sum(sym.size.reg.150[4,]<0.01)) #powers of Mardia kurt test
b5 = c(sum(sym.size.reg.150[5,]<0.1), sum(sym.size.reg.150[5,]<0.05), sum(sym.size.reg.150[5,]<0.01)) #powers of Royston test
b6 = c(sum(sym.size.reg.150[6,]<0.1), sum(sym.size.reg.150[6,]<0.05), sum(sym.size.reg.150[6,]<0.01)) #powers of MJB_M test
b7 = c(sum(sym.size.reg.150[7,]<0.1), sum(sym.size.reg.150[7,]<0.05), sum(sym.size.reg.150[7,]<0.01)) #powers of MJB_M* test
b8 = c(sum(sym.size.reg.150[8,]<0.1), sum(sym.size.reg.150[8,]<0.05), sum(sym.size.reg.150[8,]<0.01)) #powers of energy test
b9 = c(sum(sym.size.reg.150[9,]<0.1), sum(sym.size.reg.150[9,]<0.05), sum(sym.size.reg.150[9,]<0.01)) #powers of DH1 test
b10 = c(sum(sym.size.reg.150[10,]<0.1), sum(sym.size.reg.150[10,]<0.05), sum(sym.size.reg.150[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.150 = rbind(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
sym.size.reg.450 = replicate(1000, sym.size.reg(N = 450, df = 8, size = FALSE))
c1 = c(sum(sym.size.reg.450[1,]<0.1), sum(sym.size.reg.450[1,]<0.05), sum(sym.size.reg.450[1,]<0.01)) #powers of single test
c2 = c(sum(sym.size.reg.450[2,]<0.1), sum(sym.size.reg.450[2,]<0.05), sum(sym.size.reg.450[2,]<0.01)) #powers of HZ test
c3 = c(sum(sym.size.reg.450[3,]<0.1), sum(sym.size.reg.450[3,]<0.05), sum(sym.size.reg.450[3,]<0.01)) #powers of Mardia skew test
c4 = c(sum(sym.size.reg.450[4,]<0.1), sum(sym.size.reg.450[4,]<0.05), sum(sym.size.reg.450[4,]<0.01)) #powers of Mardia kurt test
c5 = c(sum(sym.size.reg.450[5,]<0.1), sum(sym.size.reg.450[5,]<0.05), sum(sym.size.reg.450[5,]<0.01)) #powers of Royston test
c6 = c(sum(sym.size.reg.450[6,]<0.1), sum(sym.size.reg.450[6,]<0.05), sum(sym.size.reg.450[6,]<0.01)) #powers of MJB_M test
c7 = c(sum(sym.size.reg.450[7,]<0.1), sum(sym.size.reg.450[7,]<0.05), sum(sym.size.reg.450[7,]<0.01)) #powers of MJB_M* test
c8 = c(sum(sym.size.reg.450[8,]<0.1), sum(sym.size.reg.450[8,]<0.05), sum(sym.size.reg.450[8,]<0.01)) #powers of energy test
c9 = c(sum(sym.size.reg.450[9,]<0.1), sum(sym.size.reg.450[9,]<0.05), sum(sym.size.reg.450[9,]<0.01)) #powers of DH1 test
c10 = c(sum(sym.size.reg.450[10,]<0.1), sum(sym.size.reg.450[10,]<0.05), sum(sym.size.reg.450[10,]<0.01)) #powers of DH2 test
result.sym.size.reg.450 = rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10)
result.sym.size.reg = cbind(result.sym.size.reg.50, result.sym.size.reg.150, result.sym.size.reg.450)
row.names(result.sym.size.reg) = c('Jarque-Bera', "Henze-Zirkler's", 'Mardia skewness', 'Mardia kurtosis', 
                                   'Royston', 'MJB_M', 'MJB_M*', 'Energy', 'Doornik-Hansen', 'Lobato-Velasco')
latex(result.sym.size.reg, numeric.dollar = FALSE) #export to latex (in working directory)

###################################
# Symulations studies for X_t = Z_1 * t + 0.5 * Z_2 * cos(pi * t) + 0.25 * Z_3 * sin(pi * t)

# IID
sizes = sym.results() # Matrix of sizes for iid
powers.t4 = sym.results(distr = 't4') # Matrix of powers (t(4) distribution)
powers.exp = sym.results(distr = 'exponential') # Matrix of powers (exponential distribution)
powers.sn <- sym.results(distr = 'sn') # Matrix of powers (sn distribution)

latex(sizes, numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.t4, numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.exp, numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.sn, numeric.dollar = FALSE) # Export to latex (in working directory)

# Regression
sizes.regression = sym.results.regression()
powers.t4.regression = sym.results.regression(distr = 't4') # Matrix of powers (t(4) distribution)
powers.exp.regression = sym.results.regression(distr = 'exponential') # Matrix of powers (exponential distribution)
powers.sn.regression <- sym.results.regression(distr = 'sn') # Matrix of powers (sn distribution)

latex(sizes.regression, numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.t4.regression, numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.exp.regression, numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.sn.regression, numeric.dollar = FALSE) # Export to latex (in working directory)
