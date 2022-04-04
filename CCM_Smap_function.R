# Multispatial CCM

CCM_EU<-function(x, column, niter, press=c("mid")){
  
  x<-droplevels(x[which(!is.na(x$value)),])
  
  if(length(levels(x$CountryGroup))>1 & nrow(x)>=50){
    
    long_ts<-as.data.frame(x %>% group_by(CountryGroup) %>% summarize(count=n(),ab=sum(Index),sd_i=sd(Index)))
    x<-droplevels(subset(x, CountryGroup %in% long_ts$CountryGroup[which(long_ts$count>3 & long_ts$ab>0 & long_ts$sd_i>0)]))
    long_ts2<-as.data.frame(x %>% group_by(CountryGroup) %>% summarize(count=n()))
    
    x<-data.frame(x %>% group_by(CountryGroup) %>% mutate(Index2=scale(Index),value2=scale(value)))
    #x<-data.frame(x %>% group_by(CountryGroup) %>% mutate(Index2=Index/Index[1],value2=value.std))
    #x<-droplevels(na.omit(x))
    
    ab<-ddply(x, .(CountryGroup), .fun=add_row, .parallel = F)
    Accm<-ab$Index[-nrow(ab)]
    #Accm<-ab$Abd[-nrow(ab)]
    val<-ddply(x, .(CountryGroup), .fun=add_row, .parallel = F)
    Bccm<-val$value[-nrow(val)]
    
    if(press=="short"){maxE<-2}
    if(press=="mid"){maxE<-min(c(min(long_ts2$count),4))}
    if(press=="long"){maxE<-min(c(min(long_ts2$count),6))}
    Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")
    
    for(E in 2:maxE) {
      Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
      Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
    }
    
    E_A<-which.max(na.omit(Emat[,1]))+1
    E_B<-which.max(na.omit(Emat[,2]))+1
    
    if(length(E_A)>0 & length(E_B)>0){
      signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1, predsteplist=1:10)
      signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1, predsteplist=1:10)
      
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_B)
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        
        CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_A)
        CCM_significance_test[2]<-1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test<-ccmtest(CCM_boot_B,CCM_boot_B)
        CCM_significance_test[1]<-1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        CCM_significance_test<-c(1,1)
      }
      
      result<-data.frame(a_cause_b=CCM_significance_test[1], # a = species, b= pressure
                         b_cause_a=CCM_significance_test[2])
      x<-as.data.frame(x %>% group_by(CountryGroup) %>% mutate(Index2=scale(Index),value2=scale(value)))
      #x<-as.data.frame(x %>% group_by(CountryGroup) %>% mutate(Index2=Zscore(Index),value2=Zscore(value)))
      block<-data.frame(ind=x$Index2, pressure=x$value2)
      #block<-data.frame(Country=as.character(x$CountryGroup), ind=x$Index2, pressure=x$value2)
      if(result[2]>0.05 & result[1]<=0.05){ # b does not cause a but a causes b
        block<-block[,c(2,1)] # effect from a to b
        res_smap<-smap_intb(block)
        #block[,2:3] <- block[,c("pressure","ind")]
        #res_smap<-ddply(block, .(Country),smap_intb)
        #strenght1<- (length(which(res_smap$V1>0))-length(which(res_smap$V1<0)))/length(res_smap$V1)
        strenght1<-res_smap[4]
        strenght2<-0
        min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA}
      if(result[1]>0.05 & result[2]<=0.05){ # b causes a but a does not cause b
        block<-block # effcet from b to a
        res_smap<-smap_intb(block)
        #res_smap<-ddply(block, .(Country),smap_intb)
        strenght1<-0
        #strenght2<-(length(which(res_smap$V1>0))-length(which(res_smap$V1<0)))/length(res_smap$V1)
        strenght2<-res_smap[4]
        min_str2<-res_smap[1];firstq_str2<-res_smap[2];med_str2<-res_smap[3];thirdq_str2<-res_smap[5];max_str2<-res_smap[6]}
        #min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA}
      if(result[1]<=0.05 & result[2]<=0.05){ # a causes b and b causes a
        res_smap<-smap_intb(block)
        #res_smap<-ddply(block, .(Country),smap_intb)
        #strenght2<-(length(which(res_smap$V1>0))-length(which(res_smap$V1<0)))/length(res_smap$V1) # effect from b to a 
        strenght2<-res_smap[4]
        #min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA
        min_str2<-res_smap[1];firstq_str2<-res_smap[2];med_str2<-res_smap[3];thirdq_str2<-res_smap[5];max_str2<-res_smap[6]
        block<-block[,c(2,1)]
        #block[,2:3] <- block[,c("pressure","ind")]
        #res_smap<-ddply(block, .(Country),smap_intb)
        #strenght1<-(length(which(res_smap$V1>0))-length(which(res_smap$V1<0)))/length(res_smap$V1) # effect from a to b
        strenght1<-smap_intb(block)[4]
      }
      if(result[1]>0.05 & result[2]>0.05){
        strenght1<-NA
        strenght2<-NA
        min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA
      }
    }else{
      strenght1<-strenght2<-min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-NA
      CCM_significance_test<-c(NA,NA)
    }
  }else{
    strenght1<-strenght2<-min_str2<-firstq_str2<-med_str2<-thirdq_str2<-max_str2<-E_A<-E_B<-NA
    CCM_significance_test<-c(NA,NA)
  }
  result2<-data.frame(a_cause_b=CCM_significance_test[1],
                      b_cause_a=CCM_significance_test[2],
                      strenght1, strenght2,min_str2,
                      firstq_str2,med_str2,thirdq_str2,max_str2,
                      E_A=E_A,
                      E_B=E_B)
  
  return(result2)
}

# S-map

smap_intb <- function(block){
  
  #block <- block[,2:3]
  
  # Determine the best theta
  theta.examined <- seq(0, 10, by = 0.1)
  th.test <-
    pforeach(
      i      = theta.examined,
      .c     = rbind,
      .cores = 7
    )({
      th.test0 <- block_lnlp(block, method = "s-map", tp = 1,columns = 2,
                             theta  = i, silent = T, num_neighbors = 0)
    })
  
  best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']
  
  # Perform multivariate S-map to quantify interaction strengths
  smapc.res <- block_lnlp(block, method = "s-map", tp = 1, columns=2,
                          theta  = best.th, num_neighbors = 0, silent = T,
                          save_smap_coefficients = T)
  smapc.tmp <- smapc.res[[1]]$smap_coefficients
  smapc.tmp <- as.data.frame(smapc.tmp)
  
  # Column names
  colnames(smapc.tmp) <- c("effect_of_pressure","constant") 
  smapc.tmp<-smapc.tmp[,1]
  smapc.tmp<-c(summary(na.omit(smapc.tmp)))
  #smapc.tmp<-mean(smapc.tmp[,1], na.rm=T)
  return(smapc.tmp)
}

# New functions for CCM and S-map

# Detrend data

detrend_data <- function(series){
  data <- data.frame(x=1:length(series), y=series)
  fit <- lm(y~x, data=data)
  fit_summary <- summary(fit)
  p_value <- fit_summary$coefficient[, 4]['x']
  
  if (p_value < 0.05){
    return(as.numeric(series - fit$fitted.values))
  } else {
    return(series)
  }
}

# Multispatial CCM

multisp_CCM<-function(x, niter){
  
  #x<-droplevels(df_press4[df_press4$Species=="Alauda arvensis",])
  
  x<-droplevels(x)
  
  if(length(levels(as.factor(x$country)))>1 & nrow(x)>=50){
    
    x_multi_ccm <- ddply(x, .(country), .fun=add_row, .parallel = F)
    Accm<-x_multi_ccm$Abd_std[-nrow(x_multi_ccm)]
    
    maxE<-4
    
    Emat<-matrix(nrow=maxE-1, ncol=1); colnames(Emat)<-c("A")
    
    for(E in 2:maxE) {
      Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
    }
    
    E_A<-which.max(na.omit(Emat[,1]))+1
    
    multisp_CCM_intern <- function(x_multi_ccm,Accm,pressure=c("temp","urb","hico","forest"),niter){
      
      Bccm<-x_multi_ccm[-nrow(x_multi_ccm),paste0(pressure,"_std")]
      
      maxE<-5
      
      Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")
      
      for(E in 2:maxE) {
        Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
        Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
      }
      
      E_A<-which.max(na.omit(Emat[,1]))+1
      E_B<-which.max(na.omit(Emat[,2]))+1
      
      if(length(E_A)>0 & length(E_B)>0){
        signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1, predsteplist=1:10)
        signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1, predsteplist=1:10)
        
        if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
           summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
          
          CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
          CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
          CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_B)
        }
        if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
           summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
          
          CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
          CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_A)
          CCM_significance_test[2]<-1
        }
        if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
           summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
          
          CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
          CCM_significance_test<-ccmtest(CCM_boot_B,CCM_boot_B)
          CCM_significance_test[1]<-1
        }
        if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
           summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
          CCM_significance_test<-c(1,1)
        }
      }
        
        result<-data.frame(a_cause_b=CCM_significance_test[1], # a = species, b= pressure
                           b_cause_a=CCM_significance_test[2])
        names(result) <- c(paste0("species_cause_",pressure),
                           paste0(pressure,"_cause_species"))
        return(result)
    }
    
      temp_on_sp <- multisp_CCM_intern(x_multi_ccm,Accm,"temp",niter)
      urb_on_sp <- multisp_CCM_intern(x_multi_ccm,Accm,"urb",niter)
      hico_on_sp <- multisp_CCM_intern(x_multi_ccm,Accm,"hico",niter)
      forest_on_sp <- multisp_CCM_intern(x_multi_ccm,Accm,"forest",niter)
      
      res_multiCCM <- cbind(temp_on_sp,urb_on_sp,hico_on_sp,forest_on_sp, E_A)
  
  }else{res_multiCCM <- data.frame(species_cause_temp=NA,temp_cause_species=NA,species_cause_urb=NA,
                                   urb_cause_species=NA,species_cause_hico=NA,hico_cause_species=NA,  
                                   species_cause_forest=NA,forest_cause_species=NA,E_A=NA)}
    
  return(res_multiCCM)
}

# S-map

smap_fun <- function(block, res_ccm){
  
  #block <- droplevels(df_press4[df_press4$Species=="Alauda arvensis" & df_press4$country=="Austria",c("Abd_std","temp_std","urb_std","hico_std","forest_std")])
  block <- block[,c("Abd_std","temp_std","urb_std","hico_std","forest_std")]
  sub_res_ccm <- res_ccm[1,c("temp_cause_species","urb_cause_species",
                        "hico_cause_species","forest_cause_species")]
  num_ccm_sig <- which(sub_res_ccm < 0.05)#sub_res_ccm[1,sub_res_ccm < 0.05]
  sub_res_ccm <- sub_res_ccm[num_ccm_sig]
  to_keep <- order(unlist(sub_res_ccm))
  if(res_ccm[1,"E_A"]<=length(to_keep)){
    to_keep <- to_keep[1:res_ccm[1,"E_A"]]
  }
  sub_res_ccm <- sub_res_ccm[to_keep]
  
  names(block) <- sub("_.*","",names(block))
  block_sub <- block[,c("Abd",sub("_.*","",names(sub_res_ccm)))]
  block_sub <- data.frame(apply(block_sub,2,Zscore))
  
  # Determine the best theta
  theta.examined <- seq(0, 10, by = 0.1)
  th.test <-
    pforeach(
      i      = theta.examined,
      .c     = rbind,
      .cores = 7
    )({
      th.test0 <- block_lnlp(block_sub, method = "s-map", tp = 1,columns=c(1:dim(block_sub)[2]),
                             target_column = 1,theta  = i,
                             silent = T, num_neighbors = 0)
    })
  
  best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']
  
  # Perform multivariate S-map to quantify interaction strengths
  smapc.res <- block_lnlp(block_sub, method = "s-map", tp = 1, columns=c(1:dim(block_sub)[2]),
                          target_column = 1, theta  = best.th,
                           num_neighbors = 0, silent = T,
                          save_smap_coefficients = T)
  smapc.tmp <- data.frame(smapc.res[[1]]$smap_coefficients)
  
  colnames(smapc.tmp) <- c(colnames(block_sub), "Constant")
  
  rho <- smapc.res[[1]]$stats$rho
  
  n_pred <- smapc.res[[1]]$stats$num_pred
  t <- rho*sqrt(n_pred-2)/sqrt(1-rho^2)
  if (t >= 0){
    pvalue <- 1 - pt(t, df = n_pred-2)    
  } else {
    pvalue <- pt(t, df = n_pred-2)
  }
  
  theta <- round(best.th, 2)
  
  return(list(coefficients = smapc.tmp, rho = rho, pvalue = pvalue, theta = theta, res_ccm=res_ccm))
}

smap_fun_signif <- function(data_ccm, res_ccm){
  data_ccm <- droplevels(data_ccm)
  res_ccm <- droplevels(res_ccm[res_ccm$Species==levels(as.factor(data_ccm$Species)),])
  res_smap2 <- res_ccm[,c("temp_cause_species","urb_cause_species",
                          "hico_cause_species","forest_cause_species","E_A")]
  if((res_ccm$temp_cause_species < 0.05 |
     res_ccm$urb_cause_species < 0.05 |
     res_ccm$hico_cause_species < 0.05 |
     res_ccm$forest_cause_species < 0.05) & sum(abs(diff(na.omit(data_ccm$urb_std))))!=0){
    res_smap <- smap_fun(data_ccm,res_smap2)
  }else{
    res_smap <- list(coefficients = NA, rho = NA, pvalue = NA, theta = NA, res_ccm =NA)
  }

  return(res_smap)
}
