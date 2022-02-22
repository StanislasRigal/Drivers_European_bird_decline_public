class.trajectory <- function (Y = NULL, X = NULL, dataset = NULL, interval_size = 0.5)
{
  if (is.null(Y) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
    stop("either 'dataset' or at least 'Y' and 'X' must be specified")
  }
  if (is.null(Y) == TRUE & is.null(Y) == TRUE) {
    Y <- dataset[,1]
    X <- dataset[,2]
  }else{
    if (class(Y) == "character" & class(X) == "character") {
      if (is.null(dataset) == TRUE) {
        stop("if 'Y' and 'X' are character, 'dataset' must exist")
      }else{
        Y <- dataset[, Y]
        X <- dataset[, X]
      }
    }else{
      if (!(class(Y) %in% c("numeric","integer")) == TRUE & !(class(X) %in% c("numeric","integer")) == TRUE) {stop("'Y' and 'X' must be either characters or vector but 'class' must be similar")}
    }
  }
  
  data <- data.frame(cbind(Y, X))
  data <- data[order(data$X),]                                                                      # ordering the X values
  
  if (length(X)<4){
    stop("time series length must be at least 4")
  }
  
  Y <- data$Y
  X <- data$X
  
  linear.model <- lm(Y~X)
  
  orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))                                                   # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
  # we have to perform a variable change to obtain relevant values in the X interval 
  # for first_order_coefficient, second_order_coefficient and intercept,
  # knowing that X'= alpha*X + beta 
  # and chi = eta*X'^2 + theta
  
  gammab  <-  orthogonal_polynomial$coefficients[3]
  delta  <-  orthogonal_polynomial$coefficients[2]
  epsilon  <-  orthogonal_polynomial$coefficients[1]
  
  alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[2]
  beta  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[1]
  
  eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta
  
  Y2<-Y*(max(X)-min(X))/(max(Y)-min(Y))                                                             # p2 and p3 are relevant when Y and X amplitudes are equivalent,
  # in particular when studying scaled-to-1 indices, Y and X amplitudes
  # may be very different, so we scaled the amplitudes to calculate p2 and p3 
  polynomial_orthonormal_basis<-lm(Y2~poly(X,2, raw=T))$coefficients
  
  if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){                                     # non linear case
    classification <- data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = (alpha^2)*gammab*eta,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),                    # points of interest
                                 p2 = (-polynomial_orthonormal_basis[2]+1)/(2*polynomial_orthonormal_basis[3]),
                                 p3 = (-polynomial_orthonormal_basis[2]-1)/(2*polynomial_orthonormal_basis[3]))
  }else{                                                                                            # linear case
    classification <- data.frame(first_order_coefficient = delta*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = 0,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+delta*beta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = NA,
                                 p2 = NA,
                                 p3 = NA)
  }
  
  classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared                                # retrieve the adjusted coefficient of determination
  
  # compute the derivaive at xm-delta and at xm + delta with delta being half of the input interval size
  derivative  <-  2*(classification$x_m-(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  derivative2  <-  2*(classification$x_m+(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  
  
  if(sign(derivative) != sign(derivative2)){                                                        # non consistent direction around x_m
    classification$derivative  <-  NA
    classification$intercept_derivative  <-  NA
  }else{                                                                                            # consistent direction around x_m
    classification$derivative  <-  mean(c(derivative, derivative2))
    classification$intercept_derivative  <-  (classification$second_order_coefficient*classification$x_m^2+classification$first_order_coefficient*classification$x_m+classification$intercept)-classification$x_m*classification$derivative
  }
  
  # compute the derivative of the curvature function
  classification$derivated_curvature  <-  -12*(classification$second_order_coefficient^2)*(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)*(classification$second_order_coefficient/abs(classification$second_order_coefficient))/
    ((1+(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)^2)^(2.5))
  
  if(classification$second_order_pvalue>0.05){classification$derivated_curvature <- NA}
  
  classification$direction <- NA                                                                    # classify the direction
  classification$direction[which(classification$derivative > 0)] <- "increase"
  classification$direction[which(classification$derivative < 0)] <- "decrease"
  classification$direction[which(is.na(classification$derivative))] <- "stable"
  classification$direction[which(as.numeric(classification$first_order_pvalue)>0.05 & as.numeric(classification$second_order_pvalue)>0.05)] <- "stable"
  
  classification$acceleration <- NA                                                                 # classify the acceleration
  classification$acceleration[which(classification$derivated_curvature < 0)] <- "accelerated"
  classification$acceleration[which(classification$derivated_curvature > 0)] <- "decelerated"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient < 0)] <- "concave"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient > 0)] <- "convex"
  classification$acceleration[which(is.na(classification$derivated_curvature))] <- "constant"
  
  classification$shape_class <- paste(classification$direction,                                       # give the final classification combining direction and acceleration
                                      classification$acceleration,
                                      sep="_")
  
  linear.model.summary <- summary(linear.model)                                                       # provide the linear approach results for comparison
  
  classification$linear_slope <- linear.model.summary$coefficients[2, 1]
  classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
  classification$linear_intercept <- linear.model.summary$coefficients[1, 1]
  
  classification$first_X_value <- X[1]
  classification$last_X_value <- X[length(X)]
  
  row.names(classification) <- "Y"
  
  return(classification)
}

mc_trend <- function(dataset,         # data
                     niter,           # number of MC simulations
                     ref_year=NULL,   # reference year
                     mid="mid",       # by defaultif ref_year is missing, equal to the mid year of the interval or to first year
                     correction=TRUE) # set the reference value to 100 and correct values below 0 before logtransformation
{
  
  b <- data.frame(t(rep(NA, 13)))
  attributes(b)$names <- c("second_order_coefficient",
                           "first_order_coefficient",
                           "strd_error",
                           "shape_class",
                           "intercept",
                           "p_1",
                           "p_2",
                           "p_3",
                           "second_order_pvalue",
                           "slope_p_value",
                           "slope",
                           "r.sq",
                           "ref_year")
  
  if(is.null(ref_year)){
    if(mid=="mid"){
      ref_year <- dataset$Year[round(nrow(dataset)/2)+1]
    }
    if(mid=="first"){
      ref_year <- min(dataset$Year)
    }
  }
  
  if(correction == TRUE){
    ref_value <- dataset$Index[dataset$Year == ref_year]
    if(ref_value == 1){
      dataset$Index <- 100*dataset$Index             # set reference year value to 100
      dataset$Index[which(dataset$Index <= 1)] <- 1    # set values < 1 to 1
      dataset$Index_SE[which(dataset$Index <= 1)] <- 0 # and their SE to 0
      dataset$Index_SE <- 100*dataset$Index_SE
    }
    if(ref_value!=1 & ref_value!=100){
      if(ref_value>1){
        dataset$Index <- dataset$Index/ref_value
        dataset$Index_SE <- dataset$Index_SE/ref_value
      }
      if(ref_value<1){
        stop("use 'correction = FALSE' when value of the reference year is strictly below 1")
      }
      dataset$Index <- 100*dataset$Index
      if(length(which(dataset$Index <= 1))>0){print("caution, low values corrected, if strongly decreasing or increasing trajectory, use respectively  first or last year as referential")}
      dataset$Index[which(dataset$Index <= 1)] <- 1
      dataset$Index_SE[which(dataset$Index <= 1)] <- 0
      dataset$Index_SE <- 100*dataset$Index_SE
    }
    if(ref_value == 100){
      dataset$Index[which(dataset$Index <= 1)] <- 1
      dataset$Index_SE[which(dataset$Index <= 1)] <- 0
    }
    
    dataset$sd <- dataset$Index_SE/dataset$Index                               # SE of log transformed data
    dataset$log <- log(dataset$Index)                                          # log transforme Y
    for(j in 1:nrow(dataset)){
      if(dataset$sd[j]>(dataset$log[j])){dataset$sd[j] <- (dataset$log[j])}    # set SE to value amplitude if SE > value (if not, it leads to huge values when resampling in the next loop)
    }
  }
  
  if(correction == FALSE){
    if(min(dataset$Index)<0){
      min_value <- abs(min(dataset$Index))
      dataset$Index <- dataset$Index+min_value
    }else{min_value <- 0}
    dataset$sd <- dataset$Index_SE
    dataset$log <- dataset$Index
  }
  
  for(i in 1:niter){
    
    a <- rnorm(nrow(dataset), mean=dataset$log, sd=dataset$sd)                 # simulate Y values from normal distribution (mean= original Y value, sd = original SE)
    if(correction == TRUE){
      a <- exp(a)/exp(a[which(dataset$Year == ref_year)])*100                    # set reference year value to 100 and retransform values if logtranformed
    }else{a <- a-min_value}
    
    a <- class.trajectory(a, dataset$Year)
    b[i, 1] <- a$second_order_coefficient
    b[i, 2] <- a$first_order_coefficient
    b[i, 3] <- a$strd_error
    b[i, 4] <- a$shape_class
    b[i, 5] <- a$intercept
    if(a$second_order_coefficient!=0){
      if(findInterval(a$p1,  c(min(dataset$Year), max(dataset$Year))) == 1){ # record changing point inside time series
        b[i, 6] <- a$p1}else{b[i, 6] <- NA}
      if(findInterval(a$p2,  c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 7] <- a$p2}else{b[i, 7] <- NA}
      if(findInterval(a$p3,  c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 8] <- a$p3}else{b[i, 8] <- NA}
    }else{
      b[i, 6] <- NA
      b[i, 7] <- NA
      b[i, 8] <- NA
    }
    b[i, 9] <- a$second_order_pvalue
    b[i, 10] <- a$linear_slope_pvalue
    b[i, 11] <- a$linear_slope
    b[i, 12] <- max(0,a$r.sq)
  }
  b[, 4] <- as.factor(b[, 4])
  b[, 13] <- rep(ref_year, nrow(b))
  return(b)
}

res_trend<-function(dataset,
                    niter,
                    ref_year=NULL,
                    mid="mid",
                    correction=TRUE){
  
  if(nrow(dataset)>3 & anyNA(dataset$Index_SE) == FALSE){
    
    simulated <- mc_trend(dataset, niter, ref_year, mid, correction)
    max_max<-NULL
    
    if(length(levels(simulated$shape_class))>1){
      test<-multinomial.theo.multcomp(simulated$shape_class, p = rep(1/length(levels(simulated$shape_class)),
                                                                     length(levels(simulated$shape_class))), prop=TRUE)
      if(min(test$p.value2[test$observed>test$expected])<0.05){
        max_shape <- row.names(test$p.value)[which(test$observed == max(test$observed[test$observed>test$expected]))]
        if(max_shape %in% c("increase_decelerated","stable_concave","decrease_accelerated")){
          max_max<-max_shape
          interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
          interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
          interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
          if(in.interval.lo(interval_quart1, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"decrease_accelerated"
          }
          if(in.interval.lo(interval_quart3, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"increase_decelerated"
          }
          max_shape <- c("increase_decelerated","stable_concave","decrease_accelerated")
        }
        if(max_shape %in% c("decrease_decelerated","stable_convex","increase_accelerated")){
          max_max<-max_shape
          interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
          interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
          interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
          if(in.interval.lo(interval_quart1, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"increase_accelerated"
          }
          if(in.interval.lo(interval_quart3, lo=(mean(simulated$p_1, na.rm=T)-sd(simulated$p_1, na.rm=T)), hi=(mean(simulated$p_1, na.rm=T)+sd(simulated$p_1, na.rm=T)))){
            max_max<-"decrease_decelerated"
          }
          max_shape <- c("decrease_decelerated","stable_convex","increase_accelerated")
        }
      }else{
        max_shape <- c("increase_constant","decrease_constant","stable_constant")[which.max(c(length(grep("increase",simulated$shape_class)),
                                                                                              length(grep("decrease",simulated$shape_class)),
                                                                                              length(grep("stable",simulated$shape_class))))]
      }
    }
    if(length(levels(simulated$shape_class)) == 1){max_shape <- levels(simulated$shape_class)}
    if(is.null(max_max)){max_max<-max_shape}
    
    alpha2 <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 1]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    sd_alpha2 <- sd(as.numeric(simulated[simulated$shape_class %in% max_shape, 1]))
    alpha1 <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 2]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    sd_alpha1 <- sd(as.numeric(simulated[simulated$shape_class %in% max_shape, 2]))
    inter <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 5]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    strd <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 3]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]))
    p_1 <- weighted.mean(as.numeric(simulated[simulated$shape_class==max_max, 6]),as.numeric(simulated[simulated$shape_class==max_max, 12]), na.rm=T)
    sd_p_1 <- sd(as.numeric(simulated[simulated$shape_class==max_max, 6]), na.rm=T)
    if( !is.na(p_1) && findInterval(p_1, c(min(dataset$Year), max(dataset$Year))) != 1){p_1 <- sd_p_1 <- as.numeric(NA)}
    p_2 <- weighted.mean(as.numeric(simulated[simulated$shape_class==max_max, 7]),as.numeric(simulated[simulated$shape_class==max_max, 12]), na.rm=T)
    sd_p_2 <- sd(as.numeric(simulated[simulated$shape_class==max_max, 7]), na.rm=T)
    if( !is.na(p_2) && findInterval(p_2, c(min(dataset$Year), max(dataset$Year))) != 1){p_2 <- sd_p_2 <- as.numeric(NA)}
    p_3 <- weighted.mean(as.numeric(simulated[simulated$shape_class==max_max, 8]),as.numeric(simulated[simulated$shape_class==max_max, 12]), na.rm=T)
    sd_p_3 <- sd(as.numeric(simulated[simulated$shape_class==max_max, 8]), na.rm=T)
    if( !is.na(p_3) && findInterval(p_3, c(min(dataset$Year), max(dataset$Year))) != 1){p_3 <- sd_p_3 <- as.numeric(NA)}
    second_order_pvalue <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 9]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]), na.rm=T)
    first_order_pvalue <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 10]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]), na.rm=T)
    slope <- weighted.mean(as.numeric(simulated[simulated$shape_class %in% max_shape, 11]),as.numeric(simulated[simulated$shape_class %in% max_shape, 12]), na.rm=T)
    slope_sd <- sd(as.numeric(simulated[simulated$shape_class %in% max_shape, 11]), na.rm=T)
    
    if(!is.null(max_max)){max_shape<-max_max}
    
  }else{alpha2 <- alpha1 <- sd_alpha1 <- inter <- strd <- p_1 <- sd_p_1 <- p_2 <- sd_p_2 <- p_3 <- sd_p_3 <- max_shape <- second_order_pvalue <- first_order_pvalue <- slope <- slope_sd <- NA}
  
  ref_year <- simulated[1,13]
  
  return(data.frame(alpha2, alpha1,sd_alpha1, inter, strd, p_1, sd_p_1, p_2, sd_p_2, p_3, sd_p_3,
                    second_order_pvalue, first_order_pvalue, slope, slope_sd, ref_year,max_shape = as.factor(max_shape)))
}


random_index_bis <- function(dataset, ref_year, dataset2, ref_value="Index"){
  dataset<-droplevels(dataset)
  
  a <- data.frame(matrix(NA, ncol=1, nrow=nrow(dataset)))
  base_year <- as.numeric(names(table(dataset$Year))[1])
  last_year <- as.numeric(names(table(dataset$Year))[length(table(dataset$Year))])
  
  if(anyNA(dataset$Index)){
    years <- dataset$Year[is.na(dataset$Index)]
    for(z in 1:length(years)){
      if(years[z]==base_year){
        k <- z
        while(k<=length(years) & is.na(dataset$Index[which(dataset$Year == (base_year+k))])){
          k <- k+1
        }
        year_t <- mean(dataset2$Index[which(dataset2$Year == years[z])])
        year_t_1 <- mean(dataset2$Index[which(dataset2$Year == years[k])])
        rat <- year_t/year_t_1
        dataset$Index[which(dataset$Year == years[z])] <- rat*dataset$Index[which(dataset$Year == (base_year+k))]
        dataset$Index_SE[which(dataset$Year == years[z])] <- 0
      }else{
        year_t <- mean(dataset2$Index[which(dataset2$Year == (years[z]-1))])
        year_t_1 <- mean(dataset2$Index[which(dataset2$Year == years[z])])
        rat <- year_t_1/year_t
        dataset$Index[which(dataset$Year == years[z])] <- rat*dataset$Index[which(dataset$Year == (years[z]-1))]
        dataset$Index_SE[which(dataset$Year == years[z])] <- 0
      }
    }
  }
  
  correction <- NULL
  
  if(length(dataset$Index[dataset$Year == ref_year])==0){
    ref_year<-base_year
  }
  
  if(dataset$Index[dataset$Year == last_year]>(5*dataset$Index[dataset$Year == ref_year])){
    correction <- dataset$Index[dataset$Year == last_year]
    dataset$Index <- dataset$Index/correction*100
    dataset$Index_SE <- dataset$Index_SE/correction*100
    ref_year2 <- last_year
  }
  if(dataset$Index[dataset$Year == base_year]>(5*dataset$Index[dataset$Year == ref_year])){
    correction <- dataset$Index[dataset$Year == base_year]
    dataset$Index <- dataset$Index/correction*100
    dataset$Index_SE <- dataset$Index_SE/correction*100
    ref_year2 <- base_year
  }
  
  dataset$sd <- dataset$Index_SE/dataset$Index
  dataset$log <- log(dataset$Index)
  for(j in 1:nrow(dataset)){
    if(dataset$sd[j]>(dataset$log[j])){dataset$sd[j] <- (dataset$log[j])}
  }
  a[,1] <- rnorm(nrow(dataset), mean = dataset$log, sd = dataset$sd)
  if(is.null(correction)){
    a_1 <- a[which(dataset$Year == ref_year), 1]
    a[,1] <- a[,1]-a_1+log(100)
  }else{
    a_1 <- a[which(dataset$Year == ref_year2), 1]
    a[,1] <- a[,1]-a_1+log(100)
  }
  
  if(ref_value %in% c("Abundance","Biomass")){
    a<-log(exp(a)/100*dataset$Index[dataset$Year==ref_year])
  }
  res<-data.frame(Year=dataset$Year,Index=a[,1], Abundance=dataset$Abundance)
  
  return(res)
}

mc_trend3 <-  function(dataset, ref_year,ref_value="Index"){
  dataset <-  droplevels(dataset)
  
  dataset$Index<-dataset$Index_imputed
  dataset$Index_SE<-dataset$Index_imputed_SE
  
  if(ref_value=="Abundance"){
    dataset$Index<-dataset$Abundance
    dataset$Index_SE<-dataset$Abundance_SE
  }
  
  if(ref_value=="Biomass"){
    dataset$Index<-dataset$biomass
    dataset$Index_SE<-dataset$biomass_SE
  }
  
  result <-  ddply(dataset, .(Species), .fun=random_index_bis, ref_year,subset(dataset, !(Species %in% levels(droplevels(as.factor(as.character(dataset$Species))[is.na(dataset$Index)])))),ref_value)
  result2 <- as.data.frame(result %>% group_by(Year) %>% summarize(Index_SE=sd(Index, na.rm=T), Index=mean(Index, na.rm=T)))
  
  if(ref_value %in% c("Abundance","Biomass")){
    result2 <- as.data.frame(result %>% group_by(Year) %>% summarize(Index_SE=sd(Index, na.rm=T), Index=log(sum(exp(Index), na.rm=T))))
  }
  return(result2)
}

msi_fun3 <-  function(dataset, ref_year, niter, ref_value="Index"){
  aaa <-  mc_trend3(dataset, ref_year, ref_value)
  for(i in 2:niter){
    aaa <-  rbind(aaa,mc_trend3(dataset, ref_year, ref_value))
  }
  
  aaa_finish<- as.data.frame(aaa %>% group_by(Year) %>%summarize(Index_SE=sd(Index, na.rm=T), Index=mean(Index, na.rm=T)))
  mean_msi <-  aaa_finish$Index
  sd_msi <-  aaa_finish$Index_SE
  
  aaa2 <-  data.frame(matrix(NA, ncol=length(table(dataset$Year)), nrow=niter))
  
  b <- data.frame(t(rep(NA, 12)))
  attributes(b)$names <- c("second_order_coef", "first_order_coef", "strd_error", "shape_class", "intercept", "p_1", "p_2", "p_3", "second_order_pvalue", "slope_p_value","linear_slope","r.sq")
  
  
  if(ref_value=="Abundance"){
    ab_init<-sum(as.numeric(dataset$Abundance[dataset$Year==ref_year]))
  }
  
  if(ref_value=="Biomass"){
    bio_init<-sum(as.numeric(dataset$biomass[dataset$Year==ref_year]))
  }
  
  for(i in 1:niter){
    aaa2[i,] <- rnorm(length(table(dataset$Year)), mean=mean_msi, sd=sd_msi)
    
    aaa2_mod <- exp(aaa2[i,]-aaa2[i, which(as.numeric(names(table(dataset$Year))) == ref_year)] + log(100))
    
    if(ref_value=="Abundance"){
      aaa2_mod<-aaa2_mod*ab_init/100
    }
    
    if(ref_value=="Biomass"){
      aaa2_mod<-aaa2_mod*bio_init/100
    }
    
    a <- class.trajectory(unlist(aaa2_mod), as.numeric(names(table(dataset$Year))))
    b[i, 1] <- a$second_order_coef
    b[i, 2] <- a$first_order_coef
    b[i, 3] <- a$strd_error
    b[i, 4] <- a$shape_class
    b[i, 5] <- a$intercept
    if(a$second_order_coef!=0){
      if(findInterval(a$p1, c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 6] <- a$p1}else{b[i, 6] <- NA}
      if(findInterval(a$p2, c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 7] <- a$p2}else{b[i, 7] <- NA}
      if(findInterval(a$p3, c(min(dataset$Year), max(dataset$Year))) == 1){
        b[i, 8] <- a$p3}else{b[i, 8] <- NA}
    }else{
      b[i, 6] <- NA
      b[i, 7] <- NA
      b[i, 8] <- NA
    }
    b[i, 9] <- a$second_order_pvalue
    b[i, 10] <- a$first_order_pvalue
    b[i, 11] <- a$linear_slope
    b[i, 12] <- max(0,a$r.sq)
  }
  b[, 4] <- as.factor(b[, 4])
  
  max_max<-NULL
  
  if(length(levels(b$shape_class))>1){
    test<-multinomial.theo.multcomp(b$shape_class, p = rep(1/length(levels(b$shape_class)),
                                                           length(levels(b$shape_class))), prop=TRUE)
    if(min(test$p.value2[test$observed>test$expected])<0.05){
      max_shape <- row.names(test$p.value)[which(test$observed == max(test$observed[test$observed>test$expected]))]
      if(max_shape %in% c("increase_decelerated","stable_concave","decrease_accelerated")){
        max_max<-max_shape
        interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
        interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
        interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
        if(in.interval.lo(interval_quart1, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"decrease_accelerated"
        }
        if(in.interval.lo(interval_quart3, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"increase_decelerated"
        }
        max_shape <- c("increase_decelerated","stable_concave","decrease_accelerated")
      }
      if(max_shape %in% c("decrease_decelerated","stable_convex","increase_accelerated")){
        max_max<-max_shape
        interval_mid<-(max(dataset$Year)-min(dataset$Year))/2+min(dataset$Year)
        interval_quart1<-interval_mid-(max(dataset$Year)-min(dataset$Year))*0.25
        interval_quart3<-interval_mid+(max(dataset$Year)-min(dataset$Year))*0.25
        if(in.interval.lo(interval_quart1, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"increase_accelerated"
        }
        if(in.interval.lo(interval_quart3, lo=(mean(b$p_1, na.rm=T)-sd(b$p_1, na.rm=T)), hi=(mean(b$p_1, na.rm=T)+sd(b$p_1, na.rm=T)))){
          max_max<-"decrease_decelerated"
        }
        max_shape <- c("decrease_decelerated","stable_convex","increase_accelerated")
      }
    }else{
      max_shape <- c("increase_constant","decrease_constant","stable_constant")[which.max(c(length(grep("increase",b$shape_class)),
                                                                                            length(grep("decrease",b$shape_class)),
                                                                                            length(grep("stable",b$shape_class))))]
    }
  }
  if(length(levels(b$shape_class)) == 1){max_shape <- levels(b$shape_class)}
  if(is.null(max_max)){max_max<-max_shape}
  
  alpha2 <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 1]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  sd_alpha2 <- sd(as.numeric(b[b$shape_class %in% max_shape, 1]))
  alpha1 <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 2]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  sd_alpha1 <- sd(as.numeric(b[b$shape_class %in% max_shape, 2]))
  inter <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 5]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  strd <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 3]),as.numeric(b[b$shape_class %in% max_shape, 12]))
  p_1 <- weighted.mean(as.numeric(b[b$shape_class==max_max, 6]),as.numeric(b[b$shape_class==max_max, 12]), na.rm=T)
  sd_p_1 <- sd(as.numeric(b[b$shape_class==max_max, 6]), na.rm=T)
  if( !is.na(p_1) && !in.interval.lo(p_1, lo=min(dataset$Year), hi=max(dataset$Year))){p_1 <- sd_p_1 <- as.numeric(NA)}
  p_2 <- weighted.mean(as.numeric(b[b$shape_class==max_max, 7]),as.numeric(b[b$shape_class==max_max, 12]), na.rm=T)
  sd_p_2 <- sd(as.numeric(b[b$shape_class==max_max, 7]), na.rm=T)
  if( !is.na(p_2) && !in.interval.lo(p_2, lo=min(dataset$Year), hi=max(dataset$Year))){p_2 <- sd_p_2 <- as.numeric(NA)}
  p_3 <- weighted.mean(as.numeric(b[b$shape_class==max_max, 8]),as.numeric(b[b$shape_class==max_max, 12]), na.rm=T)
  sd_p_3 <- sd(as.numeric(b[b$shape_class==max_max, 8]), na.rm=T)
  if( !is.na(p_3) && !in.interval.lo(p_3, lo=min(dataset$Year), hi=max(dataset$Year))){p_3 <- sd_p_3 <- as.numeric(NA)}
  second_order_pvalue <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 9]),as.numeric(b[b$shape_class %in% max_shape, 12]), na.rm=T)
  first_order_pvalue <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 10]),as.numeric(b[b$shape_class %in% max_shape, 12]), na.rm=T)
  slope <- weighted.mean(as.numeric(b[b$shape_class %in% max_shape, 11]),as.numeric(b[b$shape_class %in% max_shape, 12]), na.rm=T)
  slope_sd <- sd(as.numeric(b[b$shape_class %in% max_shape, 11]), na.rm=T)
  
  mean_msi_final<-apply(aaa2, 2, function(x){exp(mean(x,na.rm = T))})
  sd_msi_final<-apply(aaa2, 2, function(x){sd(x,na.rm = T)})
  sd_msi_final<-sd_msi_final*mean_msi_final
  sd_msi_final<-sd_msi_final/mean_msi_final[which(levels(as.factor(df_pop2$Year))==ref_year)]*100
  mean_msi_final<-mean_msi_final/mean_msi_final[which(levels(as.factor(df_pop2$Year))==ref_year)]*100
  
  if(ref_value=="Abundance"){
    mean_msi_final<-mean_msi_final*ab_init/100
    sd_msi_final<-sd_msi_final*ab_init/100
  }
  
  if(ref_value=="Biomass"){
    mean_msi_final<-mean_msi_final*bio_init/100
    sd_msi_final<-sd_msi_final*bio_init/100
  }
  
  if(!is.null(max_max)){max_shape<-max_max}
  return(list(msi = data.frame(mean_msi_final, sd_msi_final),
              coef = data.frame(alpha2, alpha1,sd_alpha1, inter, strd, p_1, sd_p_1, p_2, sd_p_2, p_3, sd_p_3,
                                second_order_pvalue, first_order_pvalue, slope, slope_sd, max_shape = as.factor(max_shape))))
}