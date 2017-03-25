#to calculate NNT/NNH for surv data
NNH_NNT_plot = function(surv, group, NNT = TRUE, time = NULL, col = "blue", treatment = NULL, 
                        ylim = NULL, ylab = "", xlab = "", main = "", 
                        yaxis = FALSE, xaxis = T, ylabel = NULL, new = FALSE,
                        lwd = 1, both.dir = F){
  require(survival)
  
  if(new){
    par(new = T)
  }
  
  if (TRUE){
    #insert values
    if (is.null(time)){
      time <- seq(range(surv[, 1]), 1)
    }
    if (length(unique(group)) != 2){
      warning("Please confirm that you are inputting two groups!\n")
    }
    if (is.null(treatment)){
      treatment <- unique(group)[1]
    }
    
    group <- relevel(as.factor(group), ref = treatment)
    
    #calculate NNT at each time point (HR approach or non-parametric approach)
    HR <- exp(coef(coxph(surv ~ group)))
    HR.high <- summary(coxph(surv ~ group))$conf.int[4]
    HR.low <- summary(coxph(surv ~ group))$conf.int[3]
    sign <- 1
    if (HR > 1){
      HR <- 1/HR
      if (NNT){
        cat("trying to plot NNT...\n")
        HR.high <- 1/HR.high
        HR.low <- 1/HR.low
      }else{
        sign <- -1
        cat("trying to plot NNH with negative sign...\n")
      }
    }else{
      if (NNT){
        sign <- -1
        HR.high <- 1/HR.high
        HR.low <- 1/HR.low
        cat("trying to plot NNT with negative sign...\n")
      }else{
        cat("trying to plot NNH...\n")
      }
    }
    
    NNT <- rep(NA, length(time))
    NNT.high <- rep(NA, length(time))
    NNT.low <- rep(NA, length(time))
    sf <- survfit(surv ~ group)
    for (i in 1:length(time)){
      s <- try(summary(sf, time = time[i]), silent = T)
      if (class(s) == "try-error"){
        NNT[i] <- NA
      }else{
        #s.exp <- s$surv[1]
        s.ctrl <- s$surv[2]
        #NNT[i] <- 1/(s.exp - s.ctrl)
        NNT[i] <- 1/(s.ctrl^HR - s.ctrl)
        NNT.high[i] <- 1/(s.ctrl^HR.high - s.ctrl)
        NNT.low[i] <- 1/(s.ctrl^HR.low - s.ctrl)
      }
    }
    
    #plot
    if (is.null(ylim)){
      ylim <- c(range(NNT, na.rm = T)[2], range(NNT, na.rm = T)[1])
    }
    if (is.infinite(ylim[1])){
      ylim[1] <- 2 * NNT[2]
    }
    NNT[which(is.infinite(NNT))] <- ylim[1]
    if (sign == 1 && !both.dir){
      plot(time, log(NNT), type = "l", col = col, ylim = log(ylim), yaxt = "n", xaxt = "n", 
           ylab = ylab, xlab = xlab, lwd = lwd, main = main)
      
      #add y axis
      if (yaxis){
        if (is.null(ylabel)){
          ylabel <- ceiling(seq(ylim[1], ylim[2], by = (ylim[2] - ylim[1])/10))
        }
        axis(2, at = log(ylabel), labels = ylabel)
      }
    }else{
      plot(time, sign * log(NNT), type = "l", col = col, ylim = c(-log(ylim[1]), log(ylim)[1]), yaxt = "n", xaxt = "n", 
           ylab = ylab, xlab = xlab, lwd = lwd, main = main)
      ylim <- c(-ylim[1], ylim[1])
      
      #add y axis
      if (yaxis){
        if (is.null(ylabel)){
          ylabel <- ceiling(seq(ylim[1], ylim[2], by = (ylim[2] - ylim[1])/10))
        }
        axis(2, at = c(-log(ylabel), log(ylabel)), labels = c(-ylabel, ylabel))
      }
    }
    
    
    #add x axis
    if (xaxis){
      axis(1, at = time, labels = time)
    }
    
    return (list(NNT = NNT*sign, NNT.high, NNT.low))
  }
}

