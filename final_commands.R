source('final_functions.R')

# IID ####
##  (5.1)
# percent = 85 %
sizes.51 <- sym.results(counter = 5000, param = TRUE)
latex(sizes.51[[1]], title = 'Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.51[[2]], title = 'SE.Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.t4.51 <- sym.results(distr = 't4', counter = 5000, param = TRUE)
latex(powers.t4.51[[1]], title = 'Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.t4.51[[2]], title = 'SE.Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.exp.51 <- sym.results(distr = 'exponential', counter = 5000, param = TRUE)
latex(powers.exp.51[[1]], title = 'Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.exp.51[[2]], title = 'SE.Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.sn.51 <- sym.results(distr = 'sn', counter = 5000, param = TRUE)
latex(powers.sn.51[[1]], title = 'Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.sn.51[[2]], title = 'SE.Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 75%
# sizes.51 <- sym.results(counter = 5000, param = TRUE, percent = .75)
# latex(sizes.51[[1]], title = 'Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(sizes.51[[2]], title = 'SE.Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.t4.51 <- sym.results(distr = 't4', counter = 5000, param = TRUE, percent = .75)
# latex(powers.t4.51[[1]], title = 'Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.t4.51[[2]], title = 'SE.Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.exp.51 <- sym.results(distr = 'exponential', counter = 5000, param = TRUE, percent = .75)
# latex(powers.exp.51[[1]], title = 'Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.exp.51[[2]], title = 'SE.Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.sn.51 <- sym.results(distr = 'sn', counter = 5000, param = TRUE, percent = .75)
# latex(powers.sn.51[[1]], title = 'Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.sn.51[[2]], title = 'SE.Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 95%
# sizes.51 <- sym.results(counter = 5000, param = TRUE, percent = .95)
# latex(sizes.51[[1]], title = 'Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(sizes.51[[2]], title = 'SE.Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.t4.51 <- sym.results(distr = 't4', counter = 5000, param = TRUE, percent = .95)
# latex(powers.t4.51[[1]], title = 'Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.t4.51[[2]], title = 'SE.Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.exp.51 <- sym.results(distr = 'exponential', counter = 5000, param = TRUE, percent = .95)
# latex(powers.exp.51[[1]], title = 'Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.exp.51[[2]], title = 'SE.Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.sn.51 <- sym.results(distr = 'sn', counter = 5000, param = TRUE, percent = .95)
# latex(powers.sn.51[[1]], title = 'Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.sn.51[[2]], title = 'SE.Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)

## (5.2)
# percent = 85%
sizes.52 <- sym.results(counter = 5000, param = FALSE)
latex(sizes.52[[1]], title = 'Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.52[[2]], title = 'SE.Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.t4.52 <- sym.results(distr = 't4', counter = 5000, param = FALSE)
latex(powers.t4.52[[1]], title = 'Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.t4.52[[2]], title = 'SE.Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.exp.52 <- sym.results(distr = 'exponential', counter = 5000, param = FALSE)
latex(powers.exp.52[[1]], title = 'Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.exp.52[[2]], title = 'SE.Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.sn.52 <- sym.results(distr = 'sn', counter = 5000, param = FALSE)
latex(powers.sn.52[[1]], title = 'Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.sn.52[[2]], title = 'SE.Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 75%
sizes.52 <- sym.results(counter = 5000, param = FALSE, percent = .75)
latex(sizes.52[[1]], title = 'Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.52[[2]], title = 'SE.Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.t4.52 <- sym.results(distr = 't4', counter = 5000, param = FALSE, percent = .75)
latex(powers.t4.52[[1]], title = 'Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.t4.52[[2]], title = 'SE.Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.exp.52 <- sym.results(distr = 'exponential', counter = 5000, param = FALSE, percent = .75)
latex(powers.exp.52[[1]], title = 'Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.exp.52[[2]], title = 'SE.Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.sn.52 <- sym.results(distr = 'sn', counter = 5000, param = FALSE, percent = .75)
latex(powers.sn.52[[1]], title = 'Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.sn.52[[2]], title = 'SE.Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 95%
sizes.52 <- sym.results(counter = 5000, param = FALSE, percent = .95)
latex(sizes.52[[1]], title = 'Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.52[[2]], title = 'SE.Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.t4.52 <- sym.results(distr = 't4', counter = 5000, param = FALSE, percent = .95)
latex(powers.t4.52[[1]], title = 'Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.t4.52[[2]], title = 'SE.Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.exp.52 <- sym.results(distr = 'exponential', counter = 5000, param = FALSE, percent = .95)
latex(powers.exp.52[[1]], title = 'Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.exp.52[[2]], title = 'SE.Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.sn.52 <- sym.results(distr = 'sn', counter = 5000, param = FALSE, percent = .95)
latex(powers.sn.52[[1]], title = 'Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.sn.52[[2]], title = 'SE.Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)

# Regression ####
## (5.1)
# percent = 85%
sizes.regresion.51 <- sym.regression.results(counter = 5000, param = TRUE)
latex(sizes.regresion.51[[1]], title = 'Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.regresion.51[[2]], title = 'SE.Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.t4.51 <- sym.regression.results(distr = 't4', counter = 5000, param = TRUE)
latex(powers.regression.t4.51[[1]], title = 'Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.t4.51[[2]], title = 'SE.Power.T4.51',numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.exp.51 <- sym.regression.results(distr = 'exponential', counter = 5000, param = TRUE)
latex(powers.regression.exp.51[[1]], title = 'Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.exp.51[[2]], title = 'SE.Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.sn.51 <- sym.regression.results(distr = 'sn', counter = 5000, param = TRUE)
latex(powers.regression.sn.51[[1]], title = 'Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.sn.51[[2]], title = 'SE.Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 75%
# sizes.regresion.51 <- sym.regression.results(counter = 5000, param = TRUE, percent = .75)
# latex(sizes.regresion.51[[1]], title = 'Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(sizes.regresion.51[[2]], title = 'SE.Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.t4.51 <- sym.regression.results(distr = 't4', counter = 5000, param = TRUE, percent = .75)
# latex(powers.regression.t4.51[[1]], title = 'Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.t4.51[[2]], title = 'SE.Power.T4.51',numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.exp.51 <- sym.regression.results(distr = 'exponential', counter = 5000, param = TRUE, percent = .75)
# latex(powers.regression.exp.51[[1]], title = 'Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.exp.51[[2]], title = 'SE.Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.sn.51 <- sym.regression.results(distr = 'sn', counter = 5000, param = TRUE, percent = .75)
# latex(powers.regression.sn.51[[1]], title = 'Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.sn.51[[2]], title = 'SE.Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 95%
# sizes.regresion.51 <- sym.regression.results(counter = 5000, param = TRUE, percent = .95)
# latex(sizes.regresion.51[[1]], title = 'Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(sizes.regresion.51[[2]], title = 'SE.Size.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.t4.51 <- sym.regression.results(distr = 't4', counter = 5000, param = TRUE, percent = .95)
# latex(powers.regression.t4.51[[1]], title = 'Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.t4.51[[2]], title = 'SE.Power.T4.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.exp.51 <- sym.regression.results(distr = 'exponential', counter = 5000, param = TRUE, percent = .95)
# latex(powers.regression.exp.51[[1]], title = 'Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.exp.51[[2]], title = 'SE.Power.Exp.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.sn.51 <- sym.regression.results(distr = 'sn', counter = 5000, param = TRUE, percent = .95)
# latex(powers.regression.sn.51[[1]], title = 'Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.sn.51[[2]], title = 'SE.Power.SN.51', numeric.dollar = FALSE) # Export to latex (in working directory)

## (5.2)
# percent = 85%
sizes.regresion.52 <- sym.regression.results(counter = 5000, param = FALSE)
latex(sizes.regresion.52[[1]], title = 'Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.regresion.52[[2]], title = 'SE.Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.t4.52 <- sym.regression.results(distr = 't4', counter = 5000, param = FALSE)
latex(powers.regression.t4.52[[1]], title = 'Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.t4.52[[2]], title = 'SE.Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
 
powers.regression.exp.52 <- sym.regression.results(distr = 'exponential', counter = 5000, param = FALSE)
latex(powers.regression.exp.52[[1]], title = 'Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.exp.52[[2]], title = 'SE.Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
 
powers.regression.sn.52 <- sym.regression.results(distr = 'sn', counter = 5000, param = FALSE)
latex(powers.regression.sn.52[[1]], title = 'Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.sn.52[[2]], title = 'SE.Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 75%
# sizes.regresion.52 <- sym.regression.results(counter = 5000, param = FALSE, percent = .75)
# latex(sizes.regresion.52[[1]], title = 'Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(sizes.regresion.52[[2]], title = 'SE.Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.t4.52 <- sym.regression.results(distr = 't4', counter = 5000, param = FALSE, percent = .75)
# latex(powers.regression.t4.52[[1]], title = 'Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.t4.52[[2]], title = 'SE.Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.exp.52 <- sym.regression.results(distr = 'exponential', counter = 5000, param = FALSE, percent = .75)
# latex(powers.regression.exp.52[[1]], title = 'Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.exp.52[[2]], title = 'SE.Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# 
# powers.regression.sn.52 <- sym.regression.results(distr = 'sn', counter = 5000, param = FALSE, percent = .75)
# latex(powers.regression.sn.52[[1]], title = 'Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)
# latex(powers.regression.sn.52[[2]], title = 'SE.Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)

# percent = 95%
sizes.regresion.52 <- sym.regression.results(counter = 5000, param = FALSE, percent = .95)
latex(sizes.regresion.52[[1]], title = 'Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(sizes.regresion.52[[2]], title = 'SE.Size.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.t4.52 <- sym.regression.results(distr = 't4', counter = 5000, param = FALSE, percent = .95)
latex(powers.regression.t4.52[[1]], title = 'Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.t4.52[[2]], title = 'SE.Power.T4.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.exp.52 <- sym.regression.results(distr = 'exponential', counter = 5000, param = FALSE, percent = .95)
latex(powers.regression.exp.52[[1]], title = 'Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.exp.52[[2]], title = 'SE.Power.Exp.52', numeric.dollar = FALSE) # Export to latex (in working directory)

powers.regression.sn.52 <- sym.regression.results(distr = 'sn', counter = 5000, param = FALSE, percent = .95)
latex(powers.regression.sn.52[[1]], title = 'Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)
latex(powers.regression.sn.52[[2]], title = 'SE.Power.SN.52', numeric.dollar = FALSE) # Export to latex (in working directory)

# Real data sets ####

library(refund)
library(fda.usc)
library(fdapace)

result.growthf <- real.data.normality.test(t(growth$hgtf))
result.growthm <- real.data.normality.test(t(growth$hgtm))

cairo_ps(filename = 'growth.eps', width = 6.5, height = 5.2)
par(mfrow = c(1, 2))
ts.plot(scale(growth$hgtf, scale = FALSE), gpars = list(main = 'Girls', xlab = 'Age', axes = FALSE))
axis(2)
axis(1, at = 1:31, labels = growth$age)
box()
ts.plot(scale(growth$hgtm, scale = FALSE), gpars = list(main = 'Boys', xlab = 'Age', axes = FALSE))
axis(2)
axis(1, at = 1:31, labels = growth$age)
box()
dev.off()

boxcox(`900 nm` ~ 1, data = as.data.frame (gasoline$NIR + 0.5), lambda = seq(-20, 2, by = 0.01))
result.gasoline.trans <- real.data.normality.test(1 / (gasoline$NIR + 1)) # With transform
result.gasoline.raw <- real.data.normality.test(gasoline$NIR + 1) # Without tansform

cairo_ps(filename = 'gasoline.eps', width = 6.5, height = 5.2)
par(mfrow = c(1, 2))
set.seed(1000)
sample.gas <- sample(1:nrow(gasoline$NIR), 10)
mean.gas <- colMeans(gasoline$NIR + 1)
ts.plot((t(gasoline$NIR + 1) - mean.gas)[, sample.gas], gpars = list(xlab = 'Wavelength (mm)', axes = FALSE, main = 'Raw'))
axis(2)
axis(1, at = seq(1, 401, by = 50), labels = seq(900, 1700, by = 100))
box()
mean.gas.rev <- colMeans(1 / (gasoline$NIR + 1))
ts.plot((t(1 / (gasoline$NIR + 1)) - mean.gas.rev)[, sample.gas], gpars = list(xlab = 'Wavelength (mm)', axes = FALSE, main = 'Transformed'))
axis(2)
axis(1, at = seq(1, 401, by = 50), labels = seq(900, 1700, by = 100))
box()
dev.off()

data(tecator)
result.tecator <- real.data.normality.test(tecator$absorp.fdata$data)

cairo_ps(filename = 'tecator.eps', width = 6.5, height = 5.2)
par(mfrow = c(1, 1))
ts.plot(scale(t(tecator$absorp.fdata$data), scale = FALSE), 
        gpars = list(xlab = tecator$absorp.fdata$names$xlab, 
                     ylab = tecator$absorp.fdata$names$ylab,
                     axes = FALSE))
axis(2)
axis(1, at = c(seq(1, 100, by = 10), 100), labels = seq(850, 1050, by = 20))
box()
dev.off()

data.set.healthy <- na.omit(cbind(DTI[DTI$case == 0, ]$cca)) # Healthy
data.set.MS <- na.omit(cbind(DTI[DTI$case == 1, ]$cca)) # MS
result.DTI.healthy <- real.data.normality.test(data.set.healthy) # Raw data
result.DTI.MS <- real.data.normality.test(data.set.MS) # Raw data

cairo_ps(filename = 'DTIRaw.eps', width = 6.5, height = 5.2)
par(mfrow = c(1, 2))
ts.plot(t(data.set.healthy), 
        gpars = list(main = 'Control - raw', 
                     xlab = 'Location',
                     ylab = 'Fractional anisotropy'))
ts.plot(t(data.set.MS), 
        gpars = list(main = 'MS - raw', 
                     xlab = 'Location',
                     ylab = 'Fractional anisotropy'))
dev.off()

cairo_ps(filename = 'DTIRes.eps', width = 6.5, height = 5.2)
par(mfrow = c(1, 2))
mean.healthy <- colMeans(data.set.healthy)
set.seed(1000)
sample.h <- sample(1:nrow(data.set.healthy), 10)
ts.plot((t(data.set.healthy) - mean.healthy)[, sample.h], 
        gpars = list(main = 'Control', 
                     xlab = 'Location',
                     ylab = 'Fractional anisotropy'))
mean.ms <- colMeans(data.set.MS)
set.seed(1000)
sample.ms <- sample(1:nrow(data.set.MS), 10)
ts.plot((t(data.set.MS) - mean.ms)[, sample.ms], 
        gpars = list(main = 'MS', 
                     xlab = 'Location',
                     ylab = 'Fractional anisotropy'))
dev.off()

cd4.list <- list()
cd4.list.times <- list()
cd4.list[[1]] <- na.omit(cd4[1,])
cd4.list.times[[1]] <- as.numeric(names(cd4.list[[1]]))
for (i in 2:nrow(cd4)) {
  cd4.list[[i]] <- na.omit(cd4[i, ])
  cd4.list.times[[i]] <- as.numeric(names(cd4.list[[i]]))
}
FPCAdense <- FPCA(cd4.list, cd4.list.times, optns = list(dataType = 'Sparse', FVEthreshold = 0.85))
result.cd4 <- real.data.normality.test.sparse(FPCAdense$xiEst)

result <- data.frame(round(rbind(unlist(result.growthf), unlist(result.growthm), unlist(result.gasoline.raw), 
                                 unlist(result.gasoline.trans), unlist(result.tecator), unlist(result.DTI.healthy), 
                                 unlist(result.DTI.MS), unlist(result.cd4)), 3))
rownames(result) <- c('growth - female', 'growth - male', 'gasoline - raw', 'gasoline - trans',
                      'tecator', 'DTI - control', 'DTI - MS', 'cd4')
colnames(result) <- c('Jarque-Bera', 'MJBM', 'Doornik-Hansen', 'Lobato-Velasco', 'p')
latex(result, title = 'results', numeric.dollar = FALSE) # Export to latex (in working directory)
