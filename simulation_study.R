####code for simulation study

library(ROCnReg)
library(sn)
library(nor1mix)

nrep <- 100
nkind <- 9
nd <- 200
nh <- 200
specificity <- c(0.75,0.8,0.85,0.9,0.95)
sensitivity <- c(0.75,0.8,0.85,0.9,0.95)

## function to generate data in 9 cases and record the nominal level
gen_data<-function(seed,  nd=200, nh=200,  specificity = c(0.75,0.8,0.85,0.9,0.95),
                   sensitivity =c(0.75,0.8,0.85,0.9,0.95)){
  specificity <- c(0.75,0.8,0.85,0.9,0.95)
  sensitivity <- c(0.75,0.8,0.85,0.9,0.95)
  if(seed==1){
    muh <- -0.75
    mud <- 2.5
    sigmad <- 1
    sigmah <- 1
    
    
    #threshold corresponding to specificity
    threshold <- qnorm(specificity, muh, sigmah)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1 - pnorm(threshold, mud, sigmad)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qnorm(1 - sensitivity, mud, sigmad)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pnorm(threshold_sen, muh, sigmah)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      yd[, l] <- rnorm(nd, mud, sigmad)
      yh[, l] <- rnorm(nh, muh, sigmah)
    }
    
  }else if(seed==2){
    muh <- 1.1
    mud <- 2.5
    sigmad <- 1
    sigmah <- 1
    
    
    #threshold corresponding to specificity
    threshold <- qnorm(specificity, muh, sigmah)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1 - pnorm(threshold, mud, sigmad)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qnorm(1 - sensitivity, mud, sigmad)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pnorm(threshold_sen, muh, sigmah)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      yd[, l] <- rnorm(nd, mud, sigmad)
      yh[, l] <- rnorm(nh, muh, sigmah)
    }
    
  }else if(seed==3){
    muh <- 2.2
    mud <- 2.5
    sigmad <- 1
    sigmah <- 1
    
    
    #threshold corresponding to specificity
    threshold <- qnorm(specificity, muh, sigmah)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1 - pnorm(threshold, mud, sigmad)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qnorm(1 - sensitivity, mud, sigmad)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pnorm(threshold_sen, muh, sigmah)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      yd[, l] <- rnorm(nd, mud, sigmad)
      yh[, l] <- rnorm(nh, muh, sigmah)
    }
    
  }else if(seed==4){
    alpha_h <- 3
    beta_h <- 1
    xi_d<-5
    omega_d<-2
    alpha_d<-5
    
    
    #threshold corresponding to specificity
    threshold <- qgamma(specificity, shape=alpha_h, rate=beta_h)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1 - psn(threshold, xi=xi_d, omega = omega_d,alpha = alpha_d)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qsn(1 - sensitivity, xi=xi_d, omega = omega_d,alpha = alpha_d)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pgamma(threshold_sen, shape=alpha_h, rate=beta_h)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      yd[, l] <- rsn(nd, xi=xi_d, omega = omega_d,alpha = alpha_d)
      yh[, l] <- rgamma(nh, shape=alpha_h, rate=beta_h)
    }
    
  }else if(seed==5){
    alpha_h <- 3
    beta_h <- 1
    xi_d<-3
    omega_d<-2
    alpha_d<-5
    
    
    #threshold corresponding to specificity
    threshold <- qgamma(specificity, shape=alpha_h, rate=beta_h)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1 - psn(threshold, xi=xi_d, omega = omega_d,alpha = alpha_d)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qsn(1 - sensitivity, xi=xi_d, omega = omega_d,alpha = alpha_d)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pgamma(threshold_sen, shape=alpha_h, rate=beta_h)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      yd[, l] <- rsn(nd, xi=xi_d, omega = omega_d,alpha = alpha_d)
      yh[, l] <- rgamma(nh, shape=alpha_h, rate=beta_h)
    }
    
  }else if(seed==6){
    alpha_h <- 3
    beta_h <- 1
    xi_d<-1.25
    omega_d<-2
    alpha_d<-5
    
    
    #threshold corresponding to specificity
    threshold <- qgamma(specificity, shape=alpha_h, rate=beta_h)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1 - psn(threshold, xi=xi_d, omega = omega_d,alpha = alpha_d)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qsn(1 - sensitivity, xi=xi_d, omega = omega_d,alpha = alpha_d)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pgamma(threshold_sen, shape=alpha_h, rate=beta_h)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      yd[, l] <- rsn(nd, xi=xi_d, omega = omega_d,alpha = alpha_d)
      yh[, l] <- rgamma(nh, shape=alpha_h, rate=beta_h)
    }
    
  }else if(seed==7){
    muh1 <- -2.5
    sigmah1 <- 1
    muh2 <- 0.5
    sigmah2 <- 1
    mud1 <- 2.5
    sigmad1 <- 1
    mud2 <- 5
    sigmad2 <- 1
    
    auxd <- norMix(mu = c(mud1,mud2), sigma = c(sigmad1,sigmad2), w = c(0.5,0.5))
    auxh <- norMix(mu = c(muh1,muh2), sigma = c(sigmah1,sigmah2), w = c(0.5,0.5))
    #threshold corresponding to specificity
    threshold <- qnorMix(p = specificity, obj = auxh)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1-pnorMix(q = threshold, obj = auxd)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qnorMix(p=1 - sensitivity, obj = auxd)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pnorMix(threshold_sen, obj = auxh)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      ud <- runif(nd)
      uh <- runif(nh)
      yd[,l] <- ifelse(ud < 0.5, rnorm(nd, mud1, sigmad1), rnorm(nd, mud2, sigmad2))
      yh[,l] <- ifelse(uh < 0.5, rnorm(nh, muh1, sigmah1), rnorm(nd, muh2, sigmah2))
    }
    
    
  }else if(seed==8){
    muh1 <- -1.15
    sigmah1 <- 1
    muh2 <- 1.5
    sigmah2 <- 1
    mud1 <- 1.5
    sigmad1 <- 1
    mud2 <- 3.5
    sigmad2 <- 1
    
    
    auxd <- norMix(mu = c(mud1,mud2), sigma = c(sigmad1,sigmad2), w = c(0.5,0.5))
    auxh <- norMix(mu = c(muh1,muh2), sigma = c(sigmah1,sigmah2), w = c(0.5,0.5))
    #threshold corresponding to specificity
    threshold <- qnorMix(p = specificity, obj = auxh)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1-pnorMix(q = threshold, obj = auxd)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qnorMix(p=1 - sensitivity, obj = auxd)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pnorMix(threshold_sen, obj = auxh)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      ud <- runif(nd)
      uh <- runif(nh)
      yd[,l] <- ifelse(ud < 0.5, rnorm(nd, mud1, sigmad1), rnorm(nd, mud2, sigmad2))
      yh[,l] <- ifelse(uh < 0.5, rnorm(nh, muh1, sigmah1), rnorm(nd, muh2, sigmah2))
    }
    
    
  }else if(seed==9){
    muh1 <- 0
    sigmah1 <- 1
    muh2 <- 3
    sigmah2 <- 1
    mud1 <- 0.5
    sigmad1 <- 1
    mud2 <- 3.25
    sigmad2 <- 1
    
    
    auxd <- norMix(mu = c(mud1,mud2), sigma = c(sigmad1,sigmad2), w = c(0.5,0.5))
    auxh <- norMix(mu = c(muh1,muh2), sigma = c(sigmah1,sigmah2), w = c(0.5,0.5))
    #threshold corresponding to specificity
    threshold <- qnorMix(p = specificity, obj = auxh)
    #true value for sensitivity at the threshold calculated just above
    sens_threshold <- 1-pnorMix(q = threshold, obj = auxd)
    
    ## in case we want to fix sensitivity
    
    #threshold corresponding to sensitivity
    threshold_sen <- qnorMix(p=1 - sensitivity, obj = auxd)
    #true value for specificity at the threshold calculated just above
    spec_threshold <- pnorMix(threshold_sen, obj = auxh)
    
    #simulating 100 datasets from normal distributions in each group
    nrep <- 100
    yd <- matrix(0, nrow = nd, ncol = nrep)
    yh <- matrix(0, nrow = nh, ncol = nrep)
    
    for(l in 1:nrep){
      ud <- runif(nd)
      uh <- runif(nh)
      yd[,l] <- ifelse(ud < 0.5, rnorm(nd, mud1, sigmad1), rnorm(nd, mud2, sigmad2))
      yh[,l] <- ifelse(uh < 0.5, rnorm(nh, muh1, sigmah1), rnorm(nd, muh2, sigmah2))
    }
    
    
  }
  out=list(yd=yd,yh=yh,sens_threshold=sens_threshold,spec_threshold=spec_threshold)
  return(out)
  
}


##size=50 simulation

nd <- 50
nh <- 50
ptm=proc.time()
set.seed(12563)
store=list()
for(i_seed in 1:nkind){
  #for(i_seed in 4:6){
  data=gen_data(seed = i_seed,nd=nd,nh=nh,sensitivity = sensitivity,specificity = specificity)
  yd=data$yd; yh=data$yh; sens_threshold=data$sens_threshold; spec_threshold=data$spec_threshold
  
  #now we will apply the procedure in the package that computed the sensitivity at a given FPF (1- specificity) and
  #corresponding 95% credible interval
  aux_coverage_sens <- aux_coverage_spec <- 0
  est_spec <- int_length_spec <- est_sens <- int_length_sens <-  matrix(0,nrow=5,ncol=nrep)
  
  for(l in 1:nrep){
    
    df <- data.frame("marker" =  c(yd[, l], yh[, l]), "status" = c(rep(1, nd), rep(0, nh)))
    
    #this is mainly to obtain the MCMC chains for the means variances and weights of each component for the two groups
    test_dpm <- pooledROC.dpm(marker = "marker",
                              group = "status", 
                              tag.h = 0, 
                              data = df, 
                              standardise = TRUE, 
                              p = c(0, 1), 
                              compute.WAIC = FALSE, 
                              compute.lpml = FALSE, 
                              compute.DIC = FALSE, 
                              prior.h = priorcontrol.dpm(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, 
                                                         L =10),
                              prior.d = priorcontrol.dpm(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, L =10), 
                              mcmc = mcmccontrol(nsave = 2000, nburn = 1000, nskip = 1),
                              parallel = "snow", 
                              ncpus = 6
    )
    
    #when fixing specificity (1-FPF)!
    test_dpm_fpf <- compute.threshold.pooledROC(test_dpm, criterion = "FPF", FPF = 1 - specificity,
                                                parallel = "snow",  ncpus = 6)
    est_sens[,l] <- test_dpm_fpf$TPF[,1]
    aux_coverage_sens=aux_coverage_sens+ifelse(test_dpm_fpf$TPF[,2] <= sens_threshold & test_dpm_fpf$TPF[,3] >= sens_threshold,1,0)
    int_length_sens[,l] <- test_dpm_fpf$TPF[,3] - test_dpm_fpf$TPF[,2]
    
    
    #note that everything is returned in terms of the FPF (= 1 - specificity)
    test_dpm_tpf <- compute.threshold.pooledROC(test_dpm, criterion = "TPF", TPF = sensitivity,
                                                parallel = "snow",  ncpus = 6)
    est_spec[,l] <- 1 - test_dpm_tpf$FPF[,1]
    aux_coverage_spec=aux_coverage_spec+ifelse((1 - test_dpm_tpf$FPF[,3]) <= spec_threshold & (1- test_dpm_tpf$FPF[,2]) >= spec_threshold,1,0) 
    int_length_spec[,l] <- (1 - test_dpm_tpf$FPF[,2]) - (1 - test_dpm_tpf$FPF[,3])
    
  }
  temp=list(est_sens=est_sens,aux_coverage_sens=aux_coverage_sens,int_length_sens=int_length_sens,
            est_spec=est_spec,aux_coverage_spec=aux_coverage_spec,int_length_spec=int_length_spec)
  store[[i_seed]]=temp
}

time1=proc.time()-ptm
time1


##size=100 simulation

nd <- 100
nh <- 100
ptm=proc.time()
set.seed(12563)
store2=list()
for(i_seed in 1:nkind){
  #for(i_seed in 4:6){
  data=gen_data(seed = i_seed,nd=nd,nh=nh,sensitivity = sensitivity,specificity = specificity)
  yd=data$yd; yh=data$yh; sens_threshold=data$sens_threshold; spec_threshold=data$spec_threshold
  
  #now we will apply the procedure in the package that computed the sensitivity at a given FPF (1- specificity) and
  #corresponding 95% credible interval
  aux_coverage_sens <- aux_coverage_spec <- 0
  est_spec <- int_length_spec <- est_sens <- int_length_sens <-  matrix(0,nrow=5,ncol=nrep)
  
  for(l in 1:nrep){
    
    df <- data.frame("marker" =  c(yd[, l], yh[, l]), "status" = c(rep(1, nd), rep(0, nh)))
    
    #this is mainly to obtain the MCMC chains for the means variances and weights of each component for the two groups
    test_dpm <- pooledROC.dpm(marker = "marker",
                              group = "status", 
                              tag.h = 0, 
                              data = df, 
                              standardise = TRUE, 
                              p = c(0, 1), 
                              compute.WAIC = FALSE, 
                              compute.lpml = FALSE, 
                              compute.DIC = FALSE, 
                              prior.h = priorcontrol.dpm(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, 
                                                         L =10),
                              prior.d = priorcontrol.dpm(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, L =10), 
                              mcmc = mcmccontrol(nsave = 2000, nburn = 1000, nskip = 1),
                              parallel = "snow", 
                              ncpus = 6
    )
    
    #when fixing specificity (1-FPF)!
    test_dpm_fpf <- compute.threshold.pooledROC(test_dpm, criterion = "FPF", FPF = 1 - specificity,
                                                parallel = "snow",  ncpus = 6)
    est_sens[,l] <- test_dpm_fpf$TPF[,1]
    aux_coverage_sens=aux_coverage_sens+ifelse(test_dpm_fpf$TPF[,2] <= sens_threshold & test_dpm_fpf$TPF[,3] >= sens_threshold,1,0)
    int_length_sens[,l] <- test_dpm_fpf$TPF[,3] - test_dpm_fpf$TPF[,2]
    
    
    #note that everything is returned in terms of the FPF (= 1 - specificity)
    test_dpm_tpf <- compute.threshold.pooledROC(test_dpm, criterion = "TPF", TPF = sensitivity,
                                                parallel = "snow",  ncpus = 6)
    est_spec[,l] <- 1 - test_dpm_tpf$FPF[,1]
    aux_coverage_spec=aux_coverage_spec+ifelse((1 - test_dpm_tpf$FPF[,3]) <= spec_threshold & (1- test_dpm_tpf$FPF[,2]) >= spec_threshold,1,0) 
    int_length_spec[,l] <- (1 - test_dpm_tpf$FPF[,2]) - (1 - test_dpm_tpf$FPF[,3])
    
  }
  temp=list(est_sens=est_sens,aux_coverage_sens=aux_coverage_sens,int_length_sens=int_length_sens,
            est_spec=est_spec,aux_coverage_spec=aux_coverage_spec,int_length_spec=int_length_spec)
  store2[[i_seed]]=temp
}


time2=proc.time()-ptm
time2


##size=200 simulation

nd <- 200
nh <- 200
ptm=proc.time()
set.seed(12563)
store3=list()
for(i_seed in 1:nkind){
  #for(i_seed in 4:6){
  data=gen_data(seed = i_seed,nd=nd,nh=nh,sensitivity = sensitivity,specificity = specificity)
  yd=data$yd; yh=data$yh; sens_threshold=data$sens_threshold; spec_threshold=data$spec_threshold
  
  #now we will apply the procedure in the package that computed the sensitivity at a given FPF (1- specificity) and
  #corresponding 95% credible interval
  aux_coverage_sens <- aux_coverage_spec <- 0
  est_spec <- int_length_spec <- est_sens <- int_length_sens <-  matrix(0,nrow=5,ncol=nrep)
  
  for(l in 1:nrep){
    
    df <- data.frame("marker" =  c(yd[, l], yh[, l]), "status" = c(rep(1, nd), rep(0, nh)))
    
    #this is mainly to obtain the MCMC chains for the means variances and weights of each component for the two groups
    test_dpm <- pooledROC.dpm(marker = "marker",
                              group = "status", 
                              tag.h = 0, 
                              data = df, 
                              standardise = TRUE, 
                              p = c(0, 1), 
                              compute.WAIC = FALSE, 
                              compute.lpml = FALSE, 
                              compute.DIC = FALSE, 
                              prior.h = priorcontrol.dpm(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, 
                                                         L =10),
                              prior.d = priorcontrol.dpm(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, L =10), 
                              mcmc = mcmccontrol(nsave = 2000, nburn = 1000, nskip = 1),
                              parallel = "snow", 
                              ncpus = 6
    )
    
    #when fixing specificity (1-FPF)!
    test_dpm_fpf <- compute.threshold.pooledROC(test_dpm, criterion = "FPF", FPF = 1 - specificity,
                                                parallel = "snow",  ncpus = 6)
    est_sens[,l] <- test_dpm_fpf$TPF[,1]
    aux_coverage_sens=aux_coverage_sens+ifelse(test_dpm_fpf$TPF[,2] <= sens_threshold & test_dpm_fpf$TPF[,3] >= sens_threshold,1,0)
    int_length_sens[,l] <- test_dpm_fpf$TPF[,3] - test_dpm_fpf$TPF[,2]
    
    
    #note that everything is returned in terms of the FPF (= 1 - specificity)
    test_dpm_tpf <- compute.threshold.pooledROC(test_dpm, criterion = "TPF", TPF = sensitivity,
                                                parallel = "snow",  ncpus = 6)
    est_spec[,l] <- 1 - test_dpm_tpf$FPF[,1]
    aux_coverage_spec=aux_coverage_spec+ifelse((1 - test_dpm_tpf$FPF[,3]) <= spec_threshold & (1- test_dpm_tpf$FPF[,2]) >= spec_threshold,1,0) 
    int_length_spec[,l] <- (1 - test_dpm_tpf$FPF[,2]) - (1 - test_dpm_tpf$FPF[,3])
    
  }
  temp=list(est_sens=est_sens,aux_coverage_sens=aux_coverage_sens,int_length_sens=int_length_sens,
            est_spec=est_spec,aux_coverage_spec=aux_coverage_spec,int_length_spec=int_length_spec)
  store3[[i_seed]]=temp
}

time3=proc.time()-ptm


###analyzing simulation results
##the function to create data frame of the results
dpm_analysis<-function(store,seed){
  aux_coverage_sens=store[[seed]]$aux_coverage_sens
  int_length_sens=store[[seed]]$int_length_sens
  est_sens=store[[seed]]$est_sens
  
  aux_coverage_spec=store[[seed]]$aux_coverage_spec
  int_length_spec=store[[seed]]$int_length_spec
  est_spec=store[[seed]]$est_spec
  
  coverage_sens <-  aux_coverage_sens/nrep
  mean_int_length_sens <- apply(int_length_sens,mean,MARGIN=1)
  mean_sens=apply(est_sens,mean,MARGIN=1)
  
  coverage_spec <-  aux_coverage_spec/nrep
  mean_int_length_spec <- apply(int_length_spec,mean,MARGIN=1)
  mean_spec=apply(est_spec,mean,MARGIN=1)
  OV=c(0.104,0.484,0.880,0.172,0.482,0.860,0.159,0.510,0.897)
  true_v=true_value(seed=seed)
  true_sens=true_v$sens_threshold
  true_spec=true_v$spec_threshold
  temp_sens=data.frame(case=seed,OVL=OV[seed],'Coverage of 95% CI'=coverage_sens,
                       'mean 95% CI length'=mean_int_length_sens,
                       'mean(Se)'=mean_sens,true_sens=true_sens,'Fixed specificity'=c(0.75,0.8,0.85,0.9,0.95))
  temp_spec=data.frame(case=seed,OVL=OV[seed],'Coverage of 95% CI'=coverage_spec,
                       'mean 95% CI length'=mean_int_length_spec,
                       'mean(Sp)'=mean_spec,true_spec=true_spec,'Fixed sensitivity'=c(0.75,0.8,0.85,0.9,0.95))
  return (list(sens=temp_sens,spec=temp_spec))
}

##combining the simulation results in dataframes

library(xtable)
#size=50
a1_sens=rbind(dpm_analysis(store,1)$sens,dpm_analysis(store,2)$sens,dpm_analysis(store,3)$sens)
a2_sens=rbind(dpm_analysis(store,4)$sens,dpm_analysis(store,5)$sens,dpm_analysis(store,6)$sens)
a3_sens=rbind(dpm_analysis(store,7)$sens,dpm_analysis(store,8)$sens,dpm_analysis(store,9)$sens)

a1_spec=rbind(dpm_analysis(store,1)$spec,dpm_analysis(store,2)$spec,dpm_analysis(store,3)$spec)
a2_spec=rbind(dpm_analysis(store,4)$spec,dpm_analysis(store,5)$spec,dpm_analysis(store,6)$spec)
a3_spec=rbind(dpm_analysis(store,7)$spec,dpm_analysis(store,8)$spec,dpm_analysis(store,9)$spec)

#size=100
aa1_sens=rbind(dpm_analysis(store2,1)$sens,dpm_analysis(store2,2)$sens,dpm_analysis(store2,3)$sens)
aa2_sens=rbind(dpm_analysis(store2,4)$sens,dpm_analysis(store2,5)$sens,dpm_analysis(store2,6)$sens)
aa3_sens=rbind(dpm_analysis(store2,7)$sens,dpm_analysis(store2,8)$sens,dpm_analysis(store2,9)$sens)

aa1_spec=rbind(dpm_analysis(store2,1)$spec,dpm_analysis(store2,2)$spec,dpm_analysis(store2,3)$spec)
aa2_spec=rbind(dpm_analysis(store2,4)$spec,dpm_analysis(store2,5)$spec,dpm_analysis(store2,6)$spec)
aa3_spec=rbind(dpm_analysis(store2,7)$spec,dpm_analysis(store2,8)$spec,dpm_analysis(store2,9)$spec)


#size=200
aaa1_sens=rbind(dpm_analysis(store3,1)$sens,dpm_analysis(store3,2)$sens,dpm_analysis(store3,3)$sens)
aaa2_sens=rbind(dpm_analysis(store3,4)$sens,dpm_analysis(store3,5)$sens,dpm_analysis(store3,6)$sens)
aaa3_sens=rbind(dpm_analysis(store3,7)$sens,dpm_analysis(store3,8)$sens,dpm_analysis(store3,9)$sens)

aaa1_spec=rbind(dpm_analysis(store3,1)$spec,dpm_analysis(store3,2)$spec,dpm_analysis(store3,3)$spec)
aaa2_spec=rbind(dpm_analysis(store3,4)$spec,dpm_analysis(store3,5)$spec,dpm_analysis(store3,6)$spec)
aaa3_spec=rbind(dpm_analysis(store3,7)$spec,dpm_analysis(store3,8)$spec,dpm_analysis(store3,9)$spec)
























time3