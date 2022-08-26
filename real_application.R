####code for the diabetes dataset application
library(ROCnReg)
library(sn)
library(nor1mix)
require(rjags)
require(mixtools) 
require(MCMCpack) 
library(data.table)

load("F:/sds材料/esra/dpm code/diabetes.RData")
set.seed(198)
disease_group=diabetes$marker[diabetes$status==1]
non_disease_group=diabetes$marker[diabetes$status==0]

yd=disease_group;yh=non_disease_group
nd=length(yd);nh=length(yh)
df <- data.frame("marker" =  c(yd, yh), "status" = c(rep(1, nd), rep(0, nh)))
grid.d <- seq(min(yd) - 100, max(yd) + 20, len = 300)
grid.h <- seq(min(yh) - 100, max(yh) + 20, len = 300)

## exploratory analysis
set.seed(42)
hist(yd,xlim=c(40,520),col='skyblue',border=F,ylim=c(0,0.0235),freq = FALSE,,main='',xlab="Glucose level")
hist(yh,add=T,col=scales::alpha('red',.5),border=F,freq = FALSE)
legend("topright",legend = c("Non-disesed group",  "Diseased group"), col = c("skyblue", scales::alpha('red',.5)), 
       pch=15)

boxplot(diabetes$marker~diabetes$status,names=c('Non-diseased','Diseased'),xlab='Diabetes status',ylab='Glucose level')




##DPM model for fixing sensitivity/specificity
test_dpm2 <- pooledROC.dpm(marker = "marker",
                           group = "status", 
                           tag.h = 0, 
                           data = df, 
                           standardise = TRUE, 
                           p = c(0, 1), 
                           compute.WAIC = TRUE, 
                           compute.lpml = FALSE, 
                           compute.DIC = TRUE, 
                           density = densitycontrol(compute = TRUE, grid.h = grid.h, grid.d = grid.d),
                           prior.h = priorcontrol.dpm(m0 = 0, S0 = 100, a = 0.1, b = 0.1, alpha = 1, 
                                                      L =20),
                           prior.d = priorcontrol.dpm(m0 = 0, S0 = 100, a = 0.1, b = 0.1, alpha = 1, L =20), 
                           mcmc = mcmccontrol(nsave = 2000, nburn = 1000, nskip = 1),
                           parallel = "snow", 
                           ncpus = 6
)

specificity <- c(0.75,0.8,0.85,0.9,0.95)
sensitivity <- c(0.75,0.8,0.85,0.9,0.95)
set.seed(180)
#when fixing specificity (1-FPF)!
aux_coverage_sens <- aux_coverage_spec <- 0
est_spec <- int_length_spec <- est_sens <- int_length_sens <-  rep(0,length(specificity))
test_dpm_fpf2 <- compute.threshold.pooledROC(test_dpm2, criterion = "FPF", FPF = 1 - specificity,
                                             parallel = "snow",  ncpus = 6)
est_sens <- test_dpm_fpf2$TPF[,1]
#aux_coverage_sens=aux_coverage_sens+ifelse(test_dpm_fpf$TPF <= sens_threshold & test_dpm_fpf$TPF[,3] >= sens_threshold,1,0)
sens_low=test_dpm_fpf2$TPF[,2];sens_high=test_dpm_fpf2$TPF[,3]
int_length_sens <- test_dpm_fpf2$TPF[,3] - test_dpm_fpf2$TPF[,2]


#note that everything is returned in terms of the FPF (= 1 - specificity)
test_dpm_tpf2 <- compute.threshold.pooledROC(test_dpm2, criterion = "TPF", TPF = sensitivity,
                                             parallel = "snow",  ncpus = 6)
est_spec <- 1 - test_dpm_tpf2$FPF[,1]
#aux_coverage_spec=aux_coverage_spec+ifelse((1 - test_dpm_tpf$FPF[,3]) <= spec_threshold & (1- test_dpm_tpf$FPF[,2]) >= spec_threshold,1,0) 
spec_low=(1 - test_dpm_tpf2$FPF[,3]);spec_high=(1 - test_dpm_tpf2$FPF[,2])
int_length_spec <- (1 - test_dpm_tpf2$FPF[,2]) - (1 - test_dpm_tpf2$FPF[,3])


real_ana_sens=data.frame(est_sens=est_sens,sens_low=sens_low,sens_high=sens_high,int_length_sens=int_length_sens,fixed=c(0.75,0.8,0.85,0.9,0.95))
real_ana_spec=data.frame(est_spec=est_spec,spec_low=spec_low,spec_high=spec_high,int_length_spec=int_length_spec,fixed=c(0.75,0.8,0.85,0.9,0.95))

##compare with normal model

test_dpm3 <- pooledROC.dpm(marker = "marker",
                           group = "status", 
                           tag.h = 0, 
                           data = df, 
                           standardise = TRUE, 
                           p = c(0, 1), 
                           compute.WAIC = TRUE, 
                           compute.lpml = FALSE, 
                           compute.DIC = TRUE, 
                           density = densitycontrol(compute = TRUE, grid.h = grid.h, grid.d = grid.d),
                           prior.h = priorcontrol.dpm(m0 = 0, S0 = 100, a = 0.1, b = 0.1, alpha = 1, 
                                                      L =1),
                           prior.d = priorcontrol.dpm(m0 = 0, S0 = 100, a = 0.1, b = 0.1, alpha = 1, L =1), 
                           mcmc = mcmccontrol(nsave = 2000, nburn = 1000, nskip = 1),
                           parallel = "snow", 
                           ncpus = 6
)

set.seed(180)
#when fixing specificity (1-FPF)!
aux_coverage_sens <- aux_coverage_spec <- 0
est_spec <- int_length_spec <- est_sens <- int_length_sens <-  rep(0,length(specificity))
test_dpm_fpf2 <- compute.threshold.pooledROC(test_dpm3, criterion = "FPF", FPF = 1 - specificity,
                                             parallel = "snow",  ncpus = 6)
est_sens <- test_dpm_fpf2$TPF[,1]
#aux_coverage_sens=aux_coverage_sens+ifelse(test_dpm_fpf$TPF <= sens_threshold & test_dpm_fpf$TPF[,3] >= sens_threshold,1,0)
sens_low=test_dpm_fpf2$TPF[,2];sens_high=test_dpm_fpf2$TPF[,3]
int_length_sens <- test_dpm_fpf2$TPF[,3] - test_dpm_fpf2$TPF[,2]


#note that everything is returned in terms of the FPF (= 1 - specificity)
test_dpm_tpf2 <- compute.threshold.pooledROC(test_dpm3, criterion = "TPF", TPF = sensitivity,
                                             parallel = "snow",  ncpus = 6)
est_spec <- 1 - test_dpm_tpf2$FPF[,1]
#aux_coverage_spec=aux_coverage_spec+ifelse((1 - test_dpm_tpf$FPF[,3]) <= spec_threshold & (1- test_dpm_tpf$FPF[,2]) >= spec_threshold,1,0) 
spec_low=(1 - test_dpm_tpf2$FPF[,3]);spec_high=(1 - test_dpm_tpf2$FPF[,2])
int_length_spec <- (1 - test_dpm_tpf2$FPF[,2]) - (1 - test_dpm_tpf2$FPF[,3])

real_ana_sens2=data.frame(est_sens=est_sens,sens_low=sens_low,sens_high=sens_high,int_length_sens=int_length_sens,fixed=c(0.75,0.8,0.85,0.9,0.95))
real_ana_spec2=data.frame(est_spec=est_spec,spec_low=spec_low,spec_high=spec_high,int_length_spec=int_length_spec,fixed=c(0.75,0.8,0.85,0.9,0.95))

