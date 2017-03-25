#**********************************************
#for SPRINT dataset
#**********************************************
rm(list = ls())
setwd("F:/R/nejm_challenge/sprint_pop/")
library(ggplot2)
library(survival)
library(rpart)
library(rpart.plot)

# import data
baseline <- read.csv("./data/baseline.csv")
outcome <- read.csv("./data/outcomes.csv")
safety <- read.csv("./data/safety.csv")
ID <- baseline$MASKID
m <- length(baseline$MASKID)

# construct surv info
surv <- Surv(time = outcome[, 3], event = outcome[, 2])
surv.adverse.sae <- Surv(time = safety$SAE_DAYS, event = safety$SAE_EVNT)

# include all available variables
baseline.catagory <- data.frame(INTENSIVE = as.factor(baseline$INTENSIVE),
                                FRS = as.factor(baseline$INCLUSIONFRS), 
                                SBP = baseline$SBP, 
                                DBP = baseline$DBP, 
                                NOAGENTS = as.factor(baseline$NOAGENTS == 0), 
                                SMOKE = as.factor(c("Never", "Former", "Current", NA)[baseline$SMOKE_3CAT]),
                                ASPIRIN = as.factor(baseline$ASPIRIN), 
                                SUB.CKD = as.factor(baseline$SUB_CKD),
                                SCREAT = baseline$SCREAT, 
                                RACE.BLACK = as.factor(baseline$RACE_BLACK), 
                                AGE = baseline$AGE,  
                                FEMALE = as.factor(baseline$FEMALE), 
                                SUB.CVD = as.factor(baseline$SUB_CVD), 
                                SUB.ClincalCVD = as.factor(baseline$SUB_CLINICALCVD), 
                                SUB.SubclinicalCVD = as.factor(baseline$SUB_SUBCLINICALCVD), 
                                CHR = baseline$CHR, 
                                GLUR = baseline$GLUR, 
                                HDL = baseline$HDL,  
                                TRR = baseline$TRR,  
                                UMALCR = baseline$UMALCR,  
                                BMI = baseline$BMI,  
                                STATIN = as.factor(baseline$STATIN)
)

### characteristics
table(baseline.catagory$INTENSIVE)
prop.table(table(baseline.catagory$INTENSIVE))
table(baseline.catagory$FEMALE)
prop.table(table(baseline.catagory$FEMALE))
table(baseline.catagory$AGE >= 74)
prop.table(table(baseline.catagory$AGE >= 74))
table(baseline.catagory$RACE.BLACK)
prop.table(table(baseline.catagory$RACE.BLACK))
table(baseline.catagory$SUB.ClincalCVD)
prop.table(table(baseline.catagory$SUB.ClincalCVD))
table(baseline.catagory$UMALCR >= 34)
prop.table(table(baseline.catagory$UMALCR >= 34))
table(cut(baseline.catagory$SBP, c(0, 132, 144, 1000)))
prop.table(table(cut(baseline.catagory$SBP, c(0, 132, 144, 1000))))
table(baseline.catagory$SUB.CKD)
prop.table(table(baseline.catagory$SUB.CKD))
table(baseline.catagory$SUB.SubclinicalCVD)
prop.table(table(baseline.catagory$SUB.SubclinicalCVD))
table(baseline.catagory$SMOKE)
prop.table(table(baseline.catagory$SMOKE))
mean(baseline.catagory$AGE)
sd(baseline.catagory$AGE)
mean(baseline.catagory$SBP)
sd(baseline.catagory$SBP)
mean(baseline.catagory$DBP)
sd(baseline.catagory$DBP)


### building training & cross-validation datasets
m <- dim(baseline)[1]
all <- 1:m
all.cutoff <- 400
all.bad <- which(surv[all, 1] <= all.cutoff & surv[all, 2] == 1)
all.cutoff <- 1600
all.good <- which(surv[all, 1] > all.cutoff)
length(all.bad)
length(all.good)
all.true <- c(rep("high risk", length(all.bad)), rep("low risk", length(all.good)))
n <- length(all.bad) + length(all.good)
all.select <- c(all.bad, all.good)
training <- 1:n
testing <- which(surv[, 1] > 400 & surv[, 1] <= 1600)

#setting weight for training
t <- table(all.true[training])
t
weights <- rep(1, length(training))
weights[which(all.true[training] == "high risk")] <- t[2]/t[1]

### training decision tree
fit <- rpart(all.true[training] ~ ., data = baseline.catagory[all.select[training], ], 
             control = rpart.control(cp = 0.025), weights = weights)
fit
plot(fit, margin = 0.1)
text(fit)
rpart.plot(fit, tweak = 1.3)

############### Characteristics of stratified patients by high/low risk ###############
testing <- which(surv[, 1] > 400 & surv[, 1] <= 1600)
all.predict <- predict(fit, newdata = baseline.catagory[testing, ])
high.risk <- all.predict[, 1] > 0.5
pred <- which(high.risk == T)
table(baseline.catagory$INTENSIVE[testing], high.risk)
prop.table(table(baseline.catagory$INTENSIVE[testing], high.risk), margin = 2)
chisq.test(table(baseline.catagory$INTENSIVE[testing], high.risk))

table(baseline.catagory$FEMALE[testing], high.risk)
prop.table(table(baseline.catagory$FEMALE[testing], high.risk), margin = 2)
chisq.test(table(baseline.catagory$FEMALE[testing], high.risk))
t <- table(baseline.catagory$FEMALE[testing][pred], baseline.catagory$INTENSIVE[testing][pred])
t <- table(baseline.catagory$FEMALE[testing][-pred], baseline.catagory$INTENSIVE[testing][-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

# table(baseline.catagory$AGE >= 74)
# prop.table(table(baseline.catagory$AGE >= 74))

table(baseline.catagory$RACE.BLACK[testing], high.risk)
prop.table(table(baseline.catagory$RACE.BLACK[testing], high.risk), margin = 2)
chisq.test(table(baseline.catagory$RACE.BLACK[testing], high.risk))
t <- table(baseline.catagory$RACE.BLACK[testing][pred], baseline.catagory$INTENSIVE[testing][pred])
t <- table(baseline.catagory$RACE.BLACK[testing][-pred], baseline.catagory$INTENSIVE[testing][-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

# table(baseline.catagory$SUB.ClincalCVD)
# prop.table(table(baseline.catagory$SUB.ClincalCVD))

# table(baseline.catagory$UMALCR >= 34)
# prop.table(table(baseline.catagory$UMALCR >= 34))

table(cut(baseline.catagory$SBP[testing], c(0, 132, 144, 1000)), high.risk)
prop.table(table(cut(baseline.catagory$SBP[testing], c(0, 132, 144, 1000)), high.risk), margin = 2)
chisq.test(table(baseline.catagory$RACE.BLACK[testing], high.risk))
t <- table(cut(baseline.catagory$SBP[testing][pred], c(0, 132, 144, 1000)), baseline.catagory$INTENSIVE[testing][pred])
t <- table(cut(baseline.catagory$SBP[testing][-pred], c(0, 132, 144, 1000)), baseline.catagory$INTENSIVE[testing][-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

table(baseline.catagory$SUB.CKD[testing], high.risk)
prop.table(table(baseline.catagory$SUB.CKD[testing], high.risk), margin = 2)
chisq.test(table(baseline.catagory$SUB.CKD[testing], high.risk))
t <- table(baseline.catagory$SUB.CKD[testing][pred], baseline.catagory$INTENSIVE[testing][pred])
t <- table(baseline.catagory$SUB.CKD[testing][-pred], baseline.catagory$INTENSIVE[testing][-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

table(baseline.catagory$SUB.SubclinicalCVD[testing], high.risk)
prop.table(table(baseline.catagory$SUB.SubclinicalCVD[testing], high.risk), margin = 2)
chisq.test(table(baseline.catagory$SUB.SubclinicalCVD[testing], high.risk))
t <- table(baseline.catagory$SUB.SubclinicalCVD[testing][pred], baseline.catagory$INTENSIVE[testing][pred])
t <- table(baseline.catagory$SUB.SubclinicalCVD[testing][-pred], baseline.catagory$INTENSIVE[testing][-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

table(baseline.catagory$SMOKE[testing], high.risk)
prop.table(table(baseline.catagory$SMOKE[testing], high.risk), margin = 2)
chisq.test(table(baseline.catagory$SMOKE[testing], high.risk))
t <- table(baseline.catagory$SMOKE[testing][pred], baseline.catagory$INTENSIVE[testing][pred])
t <- table(baseline.catagory$SMOKE[testing][-pred], baseline.catagory$INTENSIVE[testing][-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

mean(baseline.catagory$AGE[testing][pred])
sd(baseline.catagory$AGE[testing][pred])
mean(baseline.catagory$AGE[testing][-pred])
sd(baseline.catagory$AGE[testing][-pred])
t.test(baseline.catagory$AGE[testing][pred],
       baseline.catagory$AGE[testing][-pred])
mean(baseline.catagory$AGE[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)])
sd(baseline.catagory$AGE[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)])
mean(baseline.catagory$AGE[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
sd(baseline.catagory$AGE[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
mean(baseline.catagory$AGE[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)])
sd(baseline.catagory$AGE[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)])
mean(baseline.catagory$AGE[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])
sd(baseline.catagory$AGE[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])
t.test(baseline.catagory$AGE[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)],
       baseline.catagory$AGE[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
t.test(baseline.catagory$AGE[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)],
       baseline.catagory$AGE[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])

mean(baseline.catagory$SBP[testing][pred])
sd(baseline.catagory$SBP[testing][pred])
mean(baseline.catagory$SBP[testing][-pred])
sd(baseline.catagory$SBP[testing][-pred])
t.test(baseline.catagory$SBP[testing][pred],
       baseline.catagory$SBP[testing][-pred])
mean(baseline.catagory$SBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)])
sd(baseline.catagory$SBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)])
mean(baseline.catagory$SBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
sd(baseline.catagory$SBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
mean(baseline.catagory$SBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)])
sd(baseline.catagory$SBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)])
mean(baseline.catagory$SBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])
sd(baseline.catagory$SBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])
t.test(baseline.catagory$SBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)],
       baseline.catagory$SBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
t.test(baseline.catagory$SBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)],
       baseline.catagory$SBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])

mean(baseline.catagory$DBP[testing][pred])
sd(baseline.catagory$DBP[testing][pred])
mean(baseline.catagory$DBP[testing][-pred])
sd(baseline.catagory$DBP[testing][-pred])
t.test(baseline.catagory$DBP[testing][pred],
       baseline.catagory$DBP[testing][-pred])
mean(baseline.catagory$DBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)])
sd(baseline.catagory$DBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)])
mean(baseline.catagory$DBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
sd(baseline.catagory$DBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
mean(baseline.catagory$DBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)])
sd(baseline.catagory$DBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)])
mean(baseline.catagory$DBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])
sd(baseline.catagory$DBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])
t.test(baseline.catagory$DBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 1)],
       baseline.catagory$DBP[testing][pred][which(baseline.catagory$INTENSIVE[testing][pred] == 0)])
t.test(baseline.catagory$DBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 1)],
       baseline.catagory$DBP[testing][-pred][which(baseline.catagory$INTENSIVE[testing][-pred] == 0)])

############### predictive effect analysis ###########

#prognostic
all.predict <- predict(fit, newdata = baseline.catagory[testing, ])
high.risk <- all.predict[, 1] > 0.5
sf <- survfit(surv[testing, ] ~ high.risk)
summary(sf, time = 365*4)
cox <- coxph(surv[testing, ] ~ high.risk)
summary(cox)

#plot with number at risk
pdf(file = "../writing/fig.2a.prognostic_SPRINT.pdf", 5, 5)
surv.test <- surv[testing, ]
par(mar = c(10, 10, 3, 1))
sf <- survfit(surv.test ~ high.risk)
plot(sf, col = c("green", "red"), 
     ylim = c(0.85, 1.0), lwd = 1.5, cex = 0.9, 
     xlim = c(400, 1600), 
     xlab = "Time (Day)", ylab = "Event-free Probability", main = "SPRINT cross validation")
legend(x = 500, y = 0.94, 
       legend = c("low risk", "high risk"), 
       col = c("green", "red"), 
       bty = "n", lwd = 1.5, cex = 0.9)
#add patient num
to.write <- matrix(NA, nrow = 2, ncol = 7)
at <- c(400, 600, 800, 1000, 1200, 1400, 1600)
mylegend <- c("low risk", "high risk")
cex <- 0.8
for (i in 1:length(at)){
  temp <- try(summary(sf, times = at[i])$n.risk, silent = T)
  if(class(temp)=="try-error"){
    to.write[,i] <- NA
  }else{
    to.write[,i] <- summary(sf, times = at[i])$n.risk
  }
}
j <- c(1, 2)
for (i in 1:2){
  mtext(to.write[j[i],], side = 1, line = 3+i*0.8, at = at, cex = cex)
  mtext(mylegend[i], side = 1, line = 3+i*0.8, at = 100, cex = cex)
}
mtext("No. at risk", side = 1, line = 3, at = 100, cex = cex)
dev.off()


#check interaction between high risk and treatment method
all.predict <- predict(fit, newdata = baseline.catagory[testing, ])
high.risk <- all.predict[, 1] > 0.5
sf <- survfit(surv[testing, ] ~ high.risk + baseline.catagory[testing, 1])
surv.test <- surv[testing, ]
treatment.test <- baseline.catagory[testing, 1]
survdiff(surv.test[which(high.risk == T)] ~ treatment.test[which(high.risk == T)])
survdiff(surv.test[which(high.risk == F)] ~ treatment.test[which(high.risk == F)])
summary(coxph(surv.test[which(high.risk == T)] ~ treatment.test[which(high.risk == T)]))
summary(coxph(surv.test[which(high.risk == F)] ~ treatment.test[which(high.risk == F)]))
summary(coxph(surv.test ~ treatment.test * high.risk))

#plot with number at risk
pdf(file = "../writing/predictive_SPRINT_v2.pdf", 5, 5)
par(mar = c(10, 10, 3, 1))
sf <- survfit(surv.test ~ high.risk + treatment.test)
plot(sf, col = c("green", "green", "red", "red"), 
     ylim = c(0.85, 1.0), lty = c(1, 2, 1, 2), lwd = 1.5, cex = 0.9, 
     xlim = c(400, 1600), 
     xlab = "Time (Day)", ylab = "Event-free Probability")
legend(x = 500, y = 0.94, 
       legend = c("low risk, intensive", "low risk, standard", "high risk, intensive", "high risk, standard"), 
       col = c("green", "green", "red", "red"), 
       lty = c(2, 1, 2, 1), bty = "n", lwd = 1.5, cex = 0.9)
#add patient num
to.write <- matrix(NA, nrow = 4, ncol = 7)
at <- c(400, 600, 800, 1000, 1200, 1400, 1600)
mylegend <- c("low risk, intensive", "low risk, standard", "high risk, intensive", "high risk, standard")
cex <- 0.8
for (i in 1:length(at)){
  temp <- try(summary(sf, times = at[i])$n.risk, silent = T)
  if(class(temp)=="try-error"){
    to.write[,i] <- NA
  }else{
    to.write[,i] <- summary(sf, times = at[i])$n.risk
  }
}
j <- c(2, 1, 4, 3)
for (i in 1:4){
  mtext(to.write[j[i],], side = 1, line = 3+i*0.8, at = at, cex = cex)
  mtext(mylegend[i], side = 1, line = 3+i*0.8, at = 100, cex = cex)
}
mtext("No. at risk", side = 1, line = 3, at = 100, cex = cex)
dev.off()

#NNH/NNT plot
#refer to paper NNH_NNT_calculation
#NTT = 1 / {Sc(t)^h - Sc(t)}
source("./script/NNH_NNT_plot.R")
pdf(file = "../writing/NNT_SPRINT.pdf", 5, 5)
surv.test <- surv[testing, ]
treatment.test <- baseline.catagory[testing, 1]
time <- 3.26*365
NNH_NNT_plot(surv.test, group = treatment.test, 
             time = time, treatment = "1", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             xlab = "Time (Day)", ylab = "NNT", main = "SPRINT cross validation")
NNH_NNT_plot(surv.test[which(high.risk == T), ], group = treatment.test[which(high.risk == T)], 
             time = time, treatment = "1", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "red")
NNH_NNT_plot(surv.test[which(high.risk == F), ], group = treatment.test[which(high.risk == F)], 
             time = time, treatment = "1", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "green")
legend(x = 500, y = log(20), legend = c("total", "high risk", "low risk"), col = c("blue", "red", "green"), lty = 1, bty = "n", lwd = 2)
dev.off()


################## SAE analysis #################
testing <- which(surv[, 1] > 400 & surv[, 1] <= 1600)
surv.adverse.sae.test <- surv.adverse.sae[testing, ]
treatment.test <- baseline.catagory[testing, 1]
all.predict.test <- predict(fit, newdata = baseline.catagory[testing, ])
high.risk.test <- all.predict.test[, 1] > 0.5

#K-M plot
pdf(file = "../writing/SAE_SRPINT_v1.pdf", 4, 4)
sf <- survfit(surv.adverse.sae.test ~ high.risk.test + treatment.test)
summary(sf, time = 365*4)
plot(sf, col = c("green", "green", "red", "red"), lty = c(1, 2, 1, 2), lwd = 1.5, cex = 0.9, 
     xlab = "Time (Day)", ylab = "SAE-free Probability", ylim = c(0.25, 1))
legend(x = 100, y = 0.6, legend = c("low risk, intensive", "low risk, standard", "high risk, intensive", "high risk, standard"), 
       col = c("green", "green", "red", "red"), lty = c(2, 1, 2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

#cumulative hazard
pdf(file = "../writing/SAE_SRPINT_Cumhaz.pdf", 4, 4)
sf <- survfit(surv.adverse.sae.test ~ high.risk + treatment.test)
plot(sf, col = c("green", "green", "red", "red"), lty = c(1, 2, 1, 2), lwd = 1.5, cex = 0.9, 
     xlab = "Time (Day)", ylab = "Cumulative Hazard of SAE", fun = "cumhaz")
legend(x = 100, y = 1.4, legend = c("low risk, intensive", "low risk, standard", "high risk, intensive", "high risk, standard"), 
       col = c("green", "green", "red", "red"), lty = c(2, 1, 2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

survdiff(surv.adverse.sae.test[which(high.risk.test == T)] ~ treatment.test[which(high.risk.test == T)])
summary(coxph(surv.adverse.sae.test[which(high.risk.test == T)] ~ treatment.test[which(high.risk.test == T)]))
survdiff(surv.adverse.sae.test[which(high.risk.test == F)] ~ treatment.test[which(high.risk.test == F)])
summary(coxph(surv.adverse.sae.test[which(high.risk.test == F)] ~ treatment.test[which(high.risk.test == F)]))
summary(coxph(surv.adverse.sae.test ~ treatment.test * high.risk.test))

#NNH/NNT plot
#refer to paper NNH_NNT_calculation
#NTT = 1 / {Sc(t)^h - Sc(t)}
source("./script/NNH_NNT_plot.R")
pdf(file = "../writing/NNH_SPRINT.pdf", 5, 5)
### please select ###
#overall SAE
surv.adverse.sae.test <- surv.adverse.sae[testing, ]
#related SAE event
surv.adverse.sae.test <- Surv(safety$REL_SAE_DAYS, safety$REL_SAE_EVNT)[testing, ]
#hypotension
surv.adverse.sae.test <- Surv(time = safety$HYP_SAE_DAYS, event = safety$HYP_SAE_EVNT)[testing, ]
#syncope
surv.adverse.sae.test <- Surv(safety$SYN_SAE_DAYS, safety$SYN_SAE_EVNT)[testing, ]
#eletronic
surv.adverse.sae.test <- Surv(safety$ELE_SAE_DAYS, safety$ELE_SAE_EVNT)[testing, ]
#AKI
surv.adverse.sae.test <- Surv(safety$AKI_SAE_DAYS, safety$AKI_SAE_EVNT)[testing, ]
#BRA
surv.adverse.sae.test <- Surv(safety$BRA_SAE_DAYS, safety$BRA_SAE_EVNT)[testing, ]
#INJ
surv.adverse.sae.test <- Surv(safety$INJ_SAE_DAYS, safety$INJ_SAE_EVNT)[testing, ]
#end
time <- 3.26*365
treatment.test <- baseline.catagory[testing, 1]
NNH_NNT_plot(surv.adverse.sae.test, group = treatment.test, NNT = F, 
             time = time, treatment = "1", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             xlab = "Time (Day)", ylab = "NNH", main = "SPRINT cross validation")
NNH_NNT_plot(surv.adverse.sae.test[which(high.risk == T), ], group = treatment.test[which(high.risk == T)], NNT = F,
             time = time, treatment = "1", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "red")
NNH_NNT_plot(surv.adverse.sae.test[which(high.risk == F), ], group = treatment.test[which(high.risk == F)], NNT = F,
             time = time, treatment = "1", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "green")
legend(x = 1000, y = log(200), legend = c("total", "high risk", "low risk"), col = c("blue", "red", "green"), lty = 1, bty = "n", lwd = 2)
dev.off()

###plot sae ARI
pdf(file = "../../writing/SAE_SPRINT_250.pdf", 4, 4)
time <- seq(10, 1600, 10)
source("./script/NNH_NNT_plot.R")
NNH.high <- NNH_NNT_plot(surv.adverse.sae.test[which(high.risk == T), ], group = treatment.test[which(high.risk == T)], NNT = F,
             time = time, treatment = "1",
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = F, col = "red")
NNH.low <- NNH_NNT_plot(surv.adverse.sae.test[which(high.risk == F), ], group = treatment.test[which(high.risk == F)], NNT = F,
             time = time, treatment = "1",
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "green")

#plot ARI curve overtime
pdf(file = "../writing/SAE_SPRINT_ARI_Curve.pdf", 4, 4)
ARI_high <- 1/NNH.high[[1]]
ARI_low <- 1/NNH.low[[1]]
plot(time, ARI_high, col = "red", type = "l", xlab = "Time (Day)", ylab = "ARI", ylim = c(-0.01, 0.04), lwd = 2, main = "ARI in SPRINT dataset")
par(new = T)
plot(time, ARI_low, col = "green", type = "l", xlab = "Time (Day)", ylab = "ARI", ylim = c(-0.01, 0.04), lwd = 2)
abline(h = 0, col ="grey")
legend(x = 0, y = 0.04, legend = c("low risk", "high risk"), col = c("green", "red"), lty = 1, bty = "n", lwd = 2)
dev.off()

#plot ARI at different time points
# #pdf(file = "../../writing/SAE_SPRINT_250.pdf", 4, 4)
# #t <- 250, 850, 1450
# time <- c(250, 850, 1450)
# NNH.high <- NNH_NNT_plot(surv.adverse.sae.test[which(high.risk == T), ], group = treatment.test[which(high.risk == T)], NNT = F,
#              time = time, treatment = "1", 
#              yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
#              new = F, col = "red")
# NNH.low <- NNH_NNT_plot(surv.adverse.sae.test[which(high.risk == F), ], group = treatment.test[which(high.risk == F)], NNT = F,
#              time = time, treatment = "1", 
#              yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
#              new = T, col = "green")
# #plot
# #pdf(file = "../writing/SAE_SPRINT_250.pdf", 4, 4)
# ARI_high <- 1/NNH.high[[1]][1]
# ARI_high.high <- 1/NNH.high[[3]][1]
# ARI_high.low <- 1/NNH.high[[2]][1]
# ARI_low <- 1/NNH.low[[1]][1]
# ARI_low.high <- 1/NNH.low[[3]][1]
# ARI_low.low <- 1/NNH.low[[2]][1]
# barcenters <- barplot(height = c(ARI_high, ARI_low), ylim = c(-0.05, 0.07), names.arg = c("high risk", "low risk"), 
#                       ylab = "ARI", axes = T, col = "white", border = NA, xlab = "SRPINT, 250d")
# arrows(barcenters, 
#        c(ARI_high.low, ARI_low.low), 
#        barcenters,
#        c(ARI_high.high, ARI_low.high),
#        lwd = 1.5, angle = 90, code = 3, length = 0.05, col = "black")
# segments(x0 = barcenters[1]-0.2, y0 = ARI_high, x1 = barcenters[1]+0.2, y1 = ARI_high, lwd = 2)
# segments(x0 = barcenters[2]-0.2, y0 = ARI_low, x1 = barcenters[2]+0.2, y1 = ARI_low, lwd = 2)
# abline(h = 0, col ="grey")
# #calculation
# SE_high <- (ARI_high.high-ARI_high)/1.96
# SE_low <- (ARI_low.high-ARI_low)/1.96
# z <- (ARI_high-ARI_low)/sqrt(SE_high^2 + SE_low^2)
# (1 - pnorm(z))
# #dev.off()
# #pdf(file = "../writing/SAE_SPRINT_850.pdf", 4, 4)
# ARI_high <- 1/NNH.high[[1]][2]
# ARI_high.high <- 1/NNH.high[[3]][2]
# ARI_high.low <- 1/NNH.high[[2]][2]
# ARI_low <- 1/NNH.low[[1]][2]
# ARI_low.high <- 1/NNH.low[[3]][2]
# ARI_low.low <- 1/NNH.low[[2]][2]
# barcenters <- barplot(height = c(ARI_high, ARI_low), ylim = c(-0.05, 0.07), names.arg = c("high risk", "low risk"), 
#                       ylab = "ARI", axes = T, col = "white", border = NA, xlab = "SRPINT, 850d")
# arrows(barcenters, 
#        c(ARI_high.low, ARI_low.low), 
#        barcenters,
#        c(ARI_high.high, ARI_low.high),
#        lwd = 1.5, angle = 90, code = 3, length = 0.05, col = "black")
# segments(x0 = barcenters[1]-0.2, y0 = ARI_high, x1 = barcenters[1]+0.2, y1 = ARI_high, lwd = 2)
# segments(x0 = barcenters[2]-0.2, y0 = ARI_low, x1 = barcenters[2]+0.2, y1 = ARI_low, lwd = 2)
# abline(h = 0, col ="grey")
# #calculation
# SE_high <- (ARI_high.high-ARI_high)/1.96
# SE_low <- (ARI_low.high-ARI_low)/1.96
# z <- (ARI_high-ARI_low)/sqrt(SE_high^2 + SE_low^2)
# (1 - pnorm(z))
# #dev.off()
# #pdf(file = "../writing/SAE_SPRINT_1450.pdf", 4, 4)
# ARI_high <- 1/NNH.high[[1]][3]
# ARI_high.high <- 1/NNH.high[[3]][3]
# ARI_high.low <- 1/NNH.high[[2]][3]
# ARI_low <- 1/NNH.low[[1]][3]
# ARI_low.high <- 1/NNH.low[[3]][3]
# ARI_low.low <- 1/NNH.low[[2]][3]
# barcenters <- barplot(height = c(ARI_high, ARI_low), ylim = c(-0.05, 0.07), names.arg = c("high risk", "low risk"), 
#                       ylab = "ARI", axes = T, col = "white", border = NA, xlab = "SRPINT, 1450d")
# arrows(barcenters, 
#        c(ARI_high.low, ARI_low.low), 
#        barcenters,
#        c(ARI_high.high, ARI_low.high),
#        lwd = 1.5, angle = 90, code = 3, length = 0.05, col = "black")
# segments(x0 = barcenters[1]-0.2, y0 = ARI_high, x1 = barcenters[1]+0.2, y1 = ARI_high, lwd = 2)
# segments(x0 = barcenters[2]-0.2, y0 = ARI_low, x1 = barcenters[2]+0.2, y1 = ARI_low, lwd = 2)
# abline(h = 0, col ="grey")
# #calculation
# SE_high <- (ARI_high.high-ARI_high)/1.96
# SE_low <- (ARI_low.high-ARI_low)/1.96
# z <- (ARI_high-ARI_low)/sqrt(SE_high^2 + SE_low^2)
# (1 - pnorm(z))
# #dev.off()


########### plot systolic blood pressure agaist time, subgrouped by high/low risk ########
#plot ave, standard red, intensive green
pdf(file = "../writing/SBP_n_med_SPRINT.pdf", 6, 6)
ID <- baseline$MASKID
bp <- read.csv("bp.new.csv")
time <- sort(unique(bp$time))
bp$INTENSIVE <- baseline$INTENSIVE[match(bp$MASKID, baseline$MASKID)]
ID.high.risk <- baseline$MASKID[testing][which(high.risk == T)]
ID.low.risk <- baseline$MASKID[testing][which(high.risk == F)]
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 1)])[1]
  sd.intensive[i] <- sd(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 0)])[1]
  sd.standard[i] <- sd(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
row.num <- 4
mar = c(4+row.num, 8, 4, 2)
par(mar = mar)
plot(c(0, time), c(ave(baseline$SBP[which(baseline$INTENSIVE == 1)])[1], ave.intensive), 
     xlim = range(time), ylim = c(110, 150), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average SBP (mmHg)", main = "SPRINT cross validation", 
     xaxt = "n")
axis(1, at = 0:10 * 6, labels = 0:10 * 6)
par(new = T)
plot(c(0, time), c(ave(baseline$SBP[which(baseline$INTENSIVE == 0)])[1], ave.standard), 
     xlim = range(time), ylim = c(110, 150), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.high.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -6, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -12, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -12, cex = 1)
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 1)])[1]
  sd.intensive[i] <- sd(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 0)])[1]
  sd.standard[i] <- sd(bp$SBP[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
par(new = T)
plot(c(0, time), c(ave(baseline$SBP[which(baseline$INTENSIVE == 1)])[1], ave.intensive), 
     xlim = range(time), ylim = c(110, 150), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(c(0, time), c(ave(baseline$SBP[which(baseline$INTENSIVE == 0)])[1], ave.standard), 
     xlim = range(time), ylim = c(110, 150), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.low.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -12, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -12, cex = 1)
#legends
legend(x = 0, y = 150, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 25, y = 150, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

########### plot delta systolic blood pressure agaist time, subgrouped by high/low risk? ########
#z test function
z.test <- function(mu1, sd1, mu2, sd2){
  p <- pnorm((mu1-mu2)/sqrt(sd1^2 + sd2^2))
  if (p > 0.5){
    return (2 * (1-p))
  }else{
    return (2 * p)
  }
}
#plot ave, standard red, intensive green
ID <- baseline$MASKID
bp <- read.csv("bp.new.csv")
bp$SBP_delta <- rep(NA, dim(bp)[1])
#construct delta
for (i in ID){
  to.select <- which(bp$MASKID == i)
  bp$SBP_delta[to.select] <- bp$SBP[to.select] - bp$SBP[to.select][which(bp$time[to.select] == 0)]
}
time <- sort(unique(bp$time))
bp$INTENSIVE <- baseline$INTENSIVE[match(bp$MASKID, baseline$MASKID)]
ID.high.risk <- baseline$MASKID[-all.select[training]][which(high.risk == T)]
ID.low.risk <- baseline$MASKID[-all.select[training]][which(high.risk == F)]
pdf(file = "../writing/SBP_delta_n_med_SPRINT.pdf", 6, 6)
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])[1]
  sd.intensive[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])[1]
  sd.standard[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
row.num <- 4
mar = c(4+row.num, 8, 4, 2)
par(mar = mar)
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(-25, 0), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average Delta SBP (mmHg)", main = "SPRINT cross validation", 
     xaxt = "n")
axis(1, at = 0:10 * 6, labels = 0:10 * 6)
par(new = T)
plot(time, ave.standard, 
     xlim = range(time), ylim = c(-25, 0), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
write.csv(file = "../writing/temp.csv", data.frame(ave.intensive, sd.intensive/sqrt(n.intensive), ave.standard, sd.standard/sqrt(n.standard)))
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.high.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -6, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -12, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -12, cex = 1)
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])[1]
  sd.intensive[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])[1]
  sd.standard[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
write.csv(file = "../writing/temp.csv", data.frame(ave.intensive, sd.intensive/sqrt(n.intensive), ave.standard, sd.standard/sqrt(n.standard)))
par(new = T)
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(-25, 0), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(time, ave.standard, 
     xlim = range(time), ylim = c(-25, 0), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.low.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -12, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -12, cex = 1)
#legends
legend(x = 0, y = -21, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 25, y = -21, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

########### plot diastolic blood pressure agaist time, subgrouped by high/low risk? ########
#plot ave, standard red, intensive green
pdf(file = "../writing/DBP_n_med_SPRINT.pdf", 6, 6)
ID <- baseline$MASKID
bp <- read.csv("bp.new.csv")
time <- sort(unique(bp$time))
bp$INTENSIVE <- baseline$INTENSIVE[match(bp$MASKID, baseline$MASKID)]
ID.high.risk <- baseline$MASKID[-all.select[training]][which(high.risk == T)]
ID.low.risk <- baseline$MASKID[-all.select[training]][which(high.risk == F)]
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 1)])[i]
  sd.intensive[i] <- sd(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 0)])[i]
  sd.standard[i] <- sd(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
row.num <- 4
mar = c(4+row.num, 8, 4, 2)
par(mar = mar)
plot(c(0, time), c(ave(baseline$DBP[which(baseline$INTENSIVE == 1)])[1], ave.intensive), 
     xlim = range(time), ylim = c(60, 100), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average DBP (mmHg)", main = "SPRINT cross validation", 
     xaxt = "n")
axis(1, at = 0:10 * 6, labels = 0:10 * 6)
par(new = T)
plot(c(0, time), c(ave(baseline$DBP[which(baseline$INTENSIVE == 0)])[1], ave.standard), 
     xlim = range(time), ylim = c(60, 100), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.high.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -6, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -12, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -12, cex = 1)
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 1)])[1]
  sd.intensive[i] <- sd(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 0)])[1]
  sd.standard[i] <- sd(bp$DBP[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
par(new = T)
plot(c(0, time), c(ave(baseline$DBP[which(baseline$INTENSIVE == 1)])[1], ave.intensive), 
     xlim = range(time), ylim = c(60, 100), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(c(0, time), c(ave(baseline$DBP[which(baseline$INTENSIVE == 0)])[1], ave.standard), 
     xlim = range(time), ylim = c(60, 100), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.low.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -12, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -12, cex = 1)
#legends
legend(x = 0, y = 100, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 25, y = 100, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

########### plot delta diastolic blood pressure agaist time, subgrouped by high/low risk? ########
#plot ave, standard red, intensive green
ID <- baseline$MASKID
bp <- read.csv("bp.new.csv")
bp$DBP_delta <- rep(NA, dim(bp)[1])
#construct delta
for (i in ID){
  to.select <- which(bp$MASKID == i)
  bp$DBP_delta[to.select] <- bp$DBP[to.select] - bp$DBP[to.select][which(bp$time[to.select] == 0)]
}
time <- sort(unique(bp$time))
bp$INTENSIVE <- baseline$INTENSIVE[match(bp$MASKID, baseline$MASKID)]
ID.high.risk <- baseline$MASKID[-all.select[training]][which(high.risk == T)]
ID.low.risk <- baseline$MASKID[-all.select[training]][which(high.risk == F)]
pdf(file = "../writing/DBP_delta_n_med_SPRINT.pdf", 6, 6)
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])[i]
  sd.intensive[i] <- sd(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])[i]
  sd.standard[i] <- sd(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
row.num <- 4
mar = c(4+row.num, 8, 4, 2)
par(mar = mar)
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(-15, 0), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average Delta DBP (mmHg)", main = "SPRINT cross validation", 
     xaxt = "n")
axis(1, at = 0:10 * 6, labels = 0:10 * 6)
par(new = T)
plot(time, ave.standard, 
     xlim = range(time), ylim = c(-15, 0), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.high.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -6, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -12, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -12, cex = 1)
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$time == time[i] & bp$MASKID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])[1]
  sd.intensive[i] <- sd(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 1)])
  ave.standard[i] <- ave(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])[1]
  sd.standard[i] <- sd(bp$DBP_delta[to.select][which(bp$INTENSIVE[to.select] == 0)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == 1))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == 0))
}
par(new = T)
plot(c(0, time), c(ave(baseline$DBP_delta[which(baseline$INTENSIVE == 1)])[1], ave.intensive), 
     xlim = range(time), ylim = c(-15, 0), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(c(0, time), c(ave(baseline$DBP_delta[which(baseline$INTENSIVE == 0)])[1], ave.standard), 
     xlim = range(time), ylim = c(-15, 0), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard), 
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:9) * 6){
  to.select <- which(bp$time == i & bp$MASKID %in% ID.low.risk)
  med.intensive[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 1)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(bp$N_BPCLASSES[to.select][which(bp$INTENSIVE[to.select] == 0)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -12, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -12, cex = 1)
#legends
legend(x = 0, y = -12, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 25, y = -12, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

#**********************************************
#for ACCORD dataset
#**********************************************
rm(list = ls())
setwd("F:/R/nejm_challenge/ACCORD/data")
library(survival)

baseline <- read.csv("accord_key.csv", stringsAsFactors = F) #[10251, 9] #arm: Intensive BP: 1, 3; Standard BP: 2, 4
outcome <- read.csv("cvdoutcomes.csv", stringsAsFactors = F) #[10251, 34] #2: primary event; #4: primary event time
labs <- read.csv("otherlabs.csv", stringsAsFactors = F) #[118662, 11] #uacr: UMALCR
sae <- read.csv("sae.csv", stringsAsFactors = F)        #[253, 2]  #1: t <= 18m; #2: t > 18m
#only keep first sae
sae <- sae[-which(duplicated(sae[,1])), ]                      #[253, 2]
meds <- read.csv("concomitantmeds.csv")
bp <- read.csv("bloodpressure.csv")                   #[181991, 5]

###### construct data for decision tree ######
ID <- baseline$MaskID[which(baseline$arm %in% 1:2)]
ACCORDdat <- data.frame(ID = ID, 
                        GENDER = baseline$female[which(baseline$MaskID %in% ID)],
                        AGE = baseline$baseline_age[which(baseline$MaskID %in% ID)],
                        BLACK = baseline$raceclass[which(baseline$MaskID %in% ID)] == "Black",
                        subclinicalCVD = baseline$cvd_hx_baseline[which(baseline$MaskID %in% ID)],
                        PrimaryTime = outcome$fuyrs_po[which(outcome$MaskID %in% ID)],
                        PrimaryEvent = outcome$censor_po[which(outcome$MaskID %in% ID)],
                        INTENSIVE = baseline$arm[which(baseline$MaskID %in% ID)] == 1
)
baseline.labs <- labs[which(labs$Visit == "BLR"), ]
baseline.bp <- bp[which(bp$Visit == "BLR"), ]
ACCORDdat <- merge(ACCORDdat, baseline.labs, by.x = "ID", by.y = "MaskID")
ACCORDdat <- merge(ACCORDdat, baseline.bp, by.x = "ID", by.y = "MaskID", all.x = T)
colnames(ACCORDdat)

ACCORDdat <- ACCORDdat[-which(is.na(ACCORDdat$uacr)), ]
surv <- Surv(time = ACCORDdat$PrimaryTime, event = 1 - ACCORDdat$PrimaryEvent) 

### univariate analysis result
summary(coxph(surv ~ ACCORDdat$INTENSIVE))
summary(coxph(surv ~ ACCORDdat$GENDER))
summary(coxph(surv ~ ACCORDdat$AGE))
summary(coxph(surv ~ ACCORDdat$BLACK))
summary(coxph(surv ~ ACCORDdat$subclinicalCVD))
uacr <- 0.1*ACCORDdat$uacr
summary(coxph(surv ~ uacr))
summary(coxph(surv ~ ACCORDdat$sbp))
gfr <- ACCORDdat$gfr < 60
summary(coxph(surv ~ gfr))

### characteristics
table(ACCORDdat$INTENSIVE)
prop.table(table(ACCORDdat$INTENSIVE))
table(ACCORDdat$GENDER)
prop.table(table(ACCORDdat$GENDER))
table(ACCORDdat$AGE >= 74)
prop.table(table(ACCORDdat$AGE >= 74))
table(ACCORDdat$BLACK)
prop.table(table(ACCORDdat$BLACK))
table(ACCORDdat$subclinicalCVD)
prop.table(table(ACCORDdat$subclinicalCVD))
table(ACCORDdat$uacr >= 34)
prop.table(table(ACCORDdat$uacr >= 34))
table(cut(ACCORDdat$sbp, c(0, 132, 144, 1000)))
prop.table(table(cut(ACCORDdat$sbp, c(0, 132, 144, 1000))))
table(cut(ACCORDdat$dbp, c(0, 70, 1000)))
prop.table(table(cut(ACCORDdat$dbp, c(0, 70, 1000))))
table(ACCORDdat$gfr < 60)
prop.table(table(ACCORDdat$gfr < 60))
mean(ACCORDdat$AGE, na.rm = T)
sd(ACCORDdat$AGE, na.rm = T)
mean(ACCORDdat$sbp, na.rm = T)
sd(ACCORDdat$sbp, na.rm = T)
mean(ACCORDdat$dbp, na.rm = T)
sd(ACCORDdat$dbp, na.rm = T)
#z test function
z.test <- function(mu1, sd1, mu2, sd2){
  p <- pnorm((mu1-mu2)/sqrt(sd1^2 + sd2^2))
  if (p > 0.5){
    return (2 * (1-p))
  }else{
    return (2 * p)
  }
}
z.test(139.7, 15.6/sqrt(9361), 139.4, 16.1/sqrt(2258))
z.test(78.1, 11.9/sqrt(9361), 76.1, 10.3/sqrt(2258))

############### Characteristics of stratified patients by high/low risk ###############
pred <- which(ACCORDdat$AGE >= 73.5 | ACCORDdat$uacr >= 34.145 | ACCORDdat$subclinicalCVD == 1)
m <- dim(ACCORDdat)[1]
risk <- rep("low risk", m)
risk[pred] <- "high risk"
risk <- relevel(as.factor(risk), ref = "low risk")
table(ACCORDdat$INTENSIVE, risk)
prop.table(table(ACCORDdat$INTENSIVE, risk), margin = 2)
chisq.test(table(ACCORDdat$INTENSIVE, risk))

table(ACCORDdat$GENDER, risk)
prop.table(table(ACCORDdat$GENDER, risk), margin = 2)
chisq.test(table(ACCORDdat$GENDER, risk))
t <- table(ACCORDdat$GENDER[pred], ACCORDdat$INTENSIVE[pred])
t <- table(ACCORDdat$GENDER[-pred], ACCORDdat$INTENSIVE[-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

# table(baseline.catagory$AGE >= 74)
# prop.table(table(baseline.catagory$AGE >= 74))

table(ACCORDdat$BLACK, risk)
prop.table(table(ACCORDdat$BLACK, risk), margin = 2)
chisq.test(table(ACCORDdat$BLACK, risk))
t <- table(ACCORDdat$BLACK[pred], ACCORDdat$INTENSIVE[pred])
t <- table(ACCORDdat$BLACK[-pred], ACCORDdat$INTENSIVE[-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

# table(baseline.catagory$SUB.ClincalCVD)
# prop.table(table(baseline.catagory$SUB.ClincalCVD))

# table(baseline.catagory$UMALCR >= 34)
# prop.table(table(baseline.catagory$UMALCR >= 34))

table(cut(ACCORDdat$sbp, c(0, 132, 144, 1000)), risk)
prop.table(table(cut(ACCORDdat$sbp, c(0, 132, 144, 1000)), risk), margin = 2)
chisq.test(table(cut(ACCORDdat$sbp, c(0, 132, 144, 1000)), risk))
t <- table(cut(ACCORDdat$sbp[pred], c(0, 132, 144, 1000)), ACCORDdat$INTENSIVE[pred])
t <- table(cut(ACCORDdat$sbp[-pred], c(0, 132, 144, 1000)), ACCORDdat$INTENSIVE[-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

table(cut(ACCORDdat$dbp, c(0, 80, 89, 1000)), risk)
prop.table(table(cut(ACCORDdat$dbp, c(0, 80, 99, 1000)), risk), margin = 2)
chisq.test(table(cut(ACCORDdat$dbp, c(0, 80, 100, 1000)), risk))
t <- table(cut(ACCORDdat$dbp[pred], c(0, 80, 100, 1000)), ACCORDdat$INTENSIVE[pred])
t <- table(cut(ACCORDdat$dbp[-pred], c(0, 80, 100, 1000)), ACCORDdat$INTENSIVE[-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

table(gfr, risk)
prop.table(table(gfr, risk), margin = 2)
chisq.test(table(gfr, risk))
t <- table(gfr[pred], ACCORDdat$INTENSIVE[pred])
t <- table(gfr[-pred], ACCORDdat$INTENSIVE[-pred])
t
prop.table(t, margin = 2)
chisq.test(t)

mean(ACCORDdat$AGE[pred])
sd(ACCORDdat$AGE[pred])
mean(ACCORDdat$AGE[-pred])
sd(ACCORDdat$AGE[-pred])
t.test(ACCORDdat$AGE[pred],
       ACCORDdat$AGE[-pred])
mean(ACCORDdat$AGE[pred][which(ACCORDdat$INTENSIVE[pred] == 1)])
sd(ACCORDdat$AGE[pred][which(ACCORDdat$INTENSIVE[pred] == 1)])
mean(ACCORDdat$AGE[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
sd(ACCORDdat$AGE[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
mean(ACCORDdat$AGE[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)])
sd(ACCORDdat$AGE[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)])
mean(ACCORDdat$AGE[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)], na.rm = T)
sd(ACCORDdat$AGE[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)], na.rm = T)
t.test(ACCORDdat$AGE[pred][which(ACCORDdat$INTENSIVE[pred] == 1)],
       ACCORDdat$AGE[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
t.test(ACCORDdat$AGE[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)],
       ACCORDdat$AGE[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)])

mean(ACCORDdat$sbp[pred], na.rm = T)
sd(ACCORDdat$sbp[pred])
mean(ACCORDdat$sbp[-pred], na.rm = T)
sd(ACCORDdat$sbp[-pred], na.rm = T)
t.test(ACCORDdat$sbp[pred],
       ACCORDdat$sbp[-pred])
mean(ACCORDdat$sbp[pred][which(ACCORDdat$INTENSIVE[pred] == 1)])
sd(ACCORDdat$sbp[pred][which(ACCORDdat$INTENSIVE[pred] == 1)])
mean(ACCORDdat$sbp[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
sd(ACCORDdat$sbp[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
mean(ACCORDdat$sbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)])
sd(ACCORDdat$sbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)])
mean(ACCORDdat$sbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)], na.rm = T)
sd(ACCORDdat$sbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)], na.rm = T)
t.test(ACCORDdat$sbp[pred][which(ACCORDdat$INTENSIVE[pred] == 1)],
       ACCORDdat$sbp[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
t.test(ACCORDdat$sbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)],
       ACCORDdat$sbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)])

mean(ACCORDdat$dbp[pred], na.rm = T)
sd(ACCORDdat$dbp[pred])
mean(ACCORDdat$dbp[-pred], na.rm = T)
sd(ACCORDdat$dbp[-pred], na.rm = T)
t.test(ACCORDdat$dbp[pred],
       ACCORDdat$dbp[-pred])
mean(ACCORDdat$dbp[pred][which(ACCORDdat$INTENSIVE[pred] == 1)])
sd(ACCORDdat$dbp[pred][which(ACCORDdat$INTENSIVE[pred] == 1)])
mean(ACCORDdat$dbp[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
sd(ACCORDdat$dbp[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
mean(ACCORDdat$dbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)])
sd(ACCORDdat$dbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)])
mean(ACCORDdat$dbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)], na.rm = T)
sd(ACCORDdat$dbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)], na.rm = T)
t.test(ACCORDdat$dbp[pred][which(ACCORDdat$INTENSIVE[pred] == 1)],
       ACCORDdat$dbp[pred][which(ACCORDdat$INTENSIVE[pred] == 0)])
t.test(ACCORDdat$dbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 1)],
       ACCORDdat$dbp[-pred][which(ACCORDdat$INTENSIVE[-pred] == 0)])


# table(baseline.catagory$SUB.SubclinicalCVD[-all.select[training]], high.risk)
# prop.table(table(baseline.catagory$SUB.SubclinicalCVD[-all.select[training]], high.risk), margin = 2)
# chisq.test(table(baseline.catagory$SUB.SubclinicalCVD[-all.select[training]], high.risk))
# 
# table(baseline.catagory$SMOKE[-all.select[training]], high.risk)
# prop.table(table(baseline.catagory$SMOKE[-all.select[training]], high.risk), margin = 2)
# chisq.test(table(baseline.catagory$SMOKE[-all.select[training]], high.risk))


###### test decision tree performance ######
pred <- which(ACCORDdat$AGE >= 73.5 | ACCORDdat$uacr >= 34.145 | ACCORDdat$subclinicalCVD == 1)
m <- dim(ACCORDdat)[1]
risk <- rep("low risk", m)
risk[pred] <- "high risk"

#risk only
cox <- coxph(surv ~ risk)
cox
summary(cox)
#plot with number at risk
pdf(file = "../../writing/fig.2c.prognostic_ACCORD.pdf", 5, 5)
par(mar = c(10, 10, 3, 1))
surv.temp <- surv
surv.temp[, 1] <- surv[, 1]*365
sf <- survfit(surv.temp ~ risk)
plot(sf, col = c("red", "green"), 
     ylim = c(0.65, 1.0), lwd = 1.5, cex = 0.9, 
     xlim = c(0, 2600), 
     xlab = "Time (Day)", ylab = "Event-free Probability")
legend(x = 100, y = 0.85, 
       legend = c("low risk", "high risk"), 
       col = c("green", "red"), 
       bty = "n", lwd = 1.5, cex = 0.9)
#add patient num
to.write <- matrix(NA, nrow = 2, ncol = 6)
at <- c(0, 500, 1000, 1500, 2000, 2500)
mylegend <- c("low risk", "high risk")
cex <- 0.8
for (i in 1:length(at)){
  temp <- try(summary(sf, times = at[i])$n.risk, silent = T)
  if(class(temp)=="try-error"){
    to.write[,i] <- NA
  }else{
    to.write[,i] <- summary(sf, times = at[i])$n.risk
  }
}
j <- c(2, 1)
for (i in 1:2){
  mtext(to.write[j[i],], side = 1, line = 3+i*0.8, at = at, cex = cex)
  mtext(mylegend[i], side = 1, line = 3+i*0.8, at = -at[2], cex = cex)
}
mtext("No. at risk", side = 1, line = 3, at = -at[2], cex = cex)
dev.off()

#treatment only
sf <- survfit(surv ~ ACCORDdat$INTENSIVE)
plot(sf, 
     ylim = c(0.65, 1.0), lty = c(1, 2), lwd = 2, cex = 2, 
     xlim = c(0, 8), 
     xlab = "Time (Year)", ylab = "Event-free Probability")
legend(x = 1, y = 0.9, 
       legend = c("standard", "intensive"), 
       lty = c(1, 2), bty = "n", lwd = 2, cex = 1.5)
survdiff(surv ~ ACCORDdat$INTENSIVE)
cox <- coxph(surv ~ relevel(as.factor(ACCORDdat$INTENSIVE), ref = 2))
summary(cox)

#risk + treatment
sf <- survfit(surv ~ risk + ACCORDdat$INTENSIVE)
pdf(file = "../../writing/predictive_ACCORD_v2.pdf", 5, 5)
par(mar = c(10, 10, 3, 1))
surv.temp <- surv
surv.temp[, 1] <- surv[, 1]*365
sf <- survfit(surv.temp ~ risk + ACCORDdat$INTENSIVE)
plot(sf, col = c("red", "red", "green", "green"), 
     ylim = c(0.65, 1.0), lty = c(1, 2, 1, 2), lwd = 1.5, cex = 0.9, 
     xlim = c(0, 2600), 
     xlab = "Time (Day)", ylab = "Event-free Probability")
legend(x = 150, y = 0.85, 
       legend = c("low risk, intensive", "low risk, standard", "high risk, intensive", "high risk, standard"), 
       col = c("green", "green", "red", "red"), 
       lty = c(2, 1, 2, 1), bty = "n", lwd = 1.5, cex = 0.9)
#add patient num
to.write <- matrix(NA, nrow = 4, ncol = 6)
at <- c(0, 500, 1000, 1500, 2000, 2500)
mylegend <- c("low risk, intensive", "low risk, standard", "high risk, intensive", "high risk, standard")
cex <- 0.8
for (i in 1:length(at)){
  temp <- try(summary(sf, times = at[i])$n.risk, silent = T)
  if(class(temp)=="try-error"){
    to.write[,i] <- NA
  }else{
    to.write[,i] <- summary(sf, times = at[i])$n.risk
  }
}
j <- c(4, 3, 2, 1)
for (i in 1:4){
  mtext(to.write[j[i],], side = 1, line = 3+i*0.8, at = at, cex = cex)
  mtext(mylegend[i], side = 1, line = 3+i*0.8, at = -at[2], cex = cex)
}
mtext("No. at risk", side = 1, line = 3, at = -at[2], cex = cex)
dev.off()


cox <- coxph(surv ~ risk * ACCORDdat$INTENSIVE)
summary(cox)
survdiff(surv[which(risk == "high risk"), ] ~ ACCORDdat$INTENSIVE[which(risk == "high risk")])
summary(coxph(surv[which(risk == "high risk"), ] ~
                as.factor(ACCORDdat$INTENSIVE[which(risk == "high risk")])))
survdiff(surv[which(risk == "low risk"), ] ~ ACCORDdat$INTENSIVE[which(risk == "low risk")])
summary(coxph(surv[which(risk== "low risk"), ] ~
                as.factor(ACCORDdat$INTENSIVE[which(risk == "low risk")])))

#NNH/NNT plot
#refer to paper NNH_NNT_calculation
#NTT = 1 / {Sc(t)^h - Sc(t)}
source("../../sprint_pop/script/NNH_NNT_plot.R")
pdf(file = "../../writing/NNT_ACCORD.pdf", 5, 5)
time <- 3.26
par(mar = c(10, 4, 4, 4))
NNH1 <- NNH_NNT_plot(surv, group =  ACCORDdat$INTENSIVE, 
             time = time, treatment = "TRUE", 
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             xaxis = F, both.dir = T, 
             xlab = "", ylab = "NNT", main = "ACCORD cross validation")
NNH2 <- NNH_NNT_plot(surv[which(risk == "high risk"), ], group =  ACCORDdat$INTENSIVE[which(risk == "high risk")], 
             time = time, treatment = "TRUE", 
             xaxis = F, both.dir = T,
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "red")
NNH3 <- NNH_NNT_plot(surv[which(risk == "low risk"), ], group =  ACCORDdat$INTENSIVE[which(risk == "low risk")], 
             time = time, treatment = "TRUE", 
             xaxis = F, both.dir = T,
             yaxis = T, ylabel = c(20, 50, 100, 200, 500, 1000), ylim = c(1000, 20), lwd = 2,
             new = T, col = "green")
legend(x = 500, y = log(20), legend = c("total", "high risk", "low risk"), col = c("blue", "red", "green"), lty = 1, bty = "n", lwd = 2)
dev.off()

#sae
ACCORDdat$sae <- rep(0, m)
for (i in 1:dim(sae)[1]){
  j <- match(sae[i, 1], ACCORDdat$ID)
  if (!is.na(j)){
    ACCORDdat$sae[j] <- sae[i, 2]
  }
}

t <- table(strata(risk, ACCORDdat$INTENSIVE), ACCORDdat$sae > 0)
t
#ARI in high risk:
ARI_high <- t[2, 2]/sum(t[2, ]) - t[1, 2]/sum(t[1, ])
ARI_low <- t[4, 2]/sum(t[4, ]) - t[3, 2]/sum(t[3, ])
SE <- function(a, n1, c, n2){
  return (sqrt((a/n1*(1-a/n1)/n1) + c/n2*(1-c/n2)/n2))
}
SE_high <- SE(t[2, 2], sum(t[2, ]), t[1, 2], sum(t[1, ]))
SE_low <- SE(t[4, 2], sum(t[4, ]), t[3, 2], sum(t[3, ]))
z <- (ARI_high-ARI_low)/sqrt(SE_high^2 + SE_low^2)
2 * (1 - pnorm(z))
qnorm(0.975)
ARI_high - 1.96*SE_high
ARI_high + 1.96*SE_high
ARI_low - 1.96*SE_low
ARI_low + 1.96*SE_low

#NNH
#Total
tT <- table(ACCORDdat$INTENSIVE, ACCORDdat$sae > 0)
tT
ARI <- tT[2, 2]/sum(tT[2, ]) - tT[1, 2]/sum(tT[1, ])
SET <- SE(tT[2, 2], sum(tT[2, ]), tT[1, 2], sum(tT[1, ]))
1/ARI
1/(ARI + 1.96*SET)
1/(ARI - 1.96*SET)
#high-risk
1/ARI_high
1/(ARI_high + 1.96*SE_high)
1/(ARI_high - 1.96*SE_high)
#low-risk
1/ARI_low
1/(ARI_low + 1.96*SE_low)
1/(ARI_low - 1.96*SE_low)

###plot sae
# plot raw data
pdf(file = "../../writing/sae_ACCORD_rate.pdf", 8, 8)
t <- table(strata(risk, ACCORDdat$INTENSIVE), ACCORDdat$sae > 0)
t
SAE.rate <- t[, 2]/(t[, 1] + t[, 2])
SE <- function(a, n){
  return (sqrt(a/n*(1-a/n)/n))
}
SE.SAE.rate <- c(SE(t[1, 1], t[1, 1] + t[1, 2]),
                 SE(t[2, 1], t[2, 1] + t[2, 2]),
                 SE(t[3, 1], t[3, 1] + t[3, 2]),
                 SE(t[4, 1], t[4, 1] + t[4, 2]))
barcenters <- barplot(height = SAE.rate, ylim = c(0, 0.06), 
                      names.arg = c("high risk, standard", "high risk, intensive", "low risk, standard", "low risk, intensive"),
                      ylab = "SAE Rate", axes = T, 
                      col = c(rgb(1, 0, 0, 0.5), rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 1, 0, 0.5)), xlab = "ACCORD")
arrows(barcenters,
       SAE.rate - 1.96*SE.SAE.rate,
       barcenters,
       SAE.rate + 1.96*SE.SAE.rate,
       lwd = 1.5, angle = 90, code = 3, length = 0.05, col = "black")
dev.off()

########## plot systolic blood pressure agaist time, subgrouped by high/low risk? ##########
#plot ave, standard red, intensive green
pdf(file = "../../writing/SBP_n_med_ACCORD_v1.pdf", 6, 6)
par(mar = c(8, 8, 4, 2))
#load bp
bp <- read.csv("bloodpressure.csv", stringsAsFactors = F)                   #[181991, 5]
bp$Visit[which(bp$Visit == "BLR")] <- 0
bp <- bp[-which(bp$Visit == "EXIT" | is.na(bp$Visit)), ]
bp$Visit <- as.numeric(gsub("F", "", bp$Visit))
bp$INTENSIVE <- ACCORDdat$INTENSIVE[match(bp$MaskID, ACCORDdat$ID)]
ID.high.risk <- ACCORDdat[pred, "ID"]
ID.low.risk <- ACCORDdat[-pred, "ID"]
time <- sort(unique(bp$Visit))
time <- c(0:3, (1:21) * 4)
#load n_med
n_med <- read.csv("concomitantmeds.csv", stringsAsFactors = F)               #[60560, 56]
names(n_med)
unique(n_med$Visit)
summary(n_med)
n_med$Visit[which(n_med$Visit == "BLR")] <- 0
n_med$Sum <- rowSums(n_med[, c(3:5, 7:16)], na.rm = T) #select BP treatment only
n_med <- n_med[-which(n_med$Visit == "EXIT"), ]
n_med$Visit <- as.numeric(gsub("F", "", n_med$Visit))
n_med$INTENSIVE <- ACCORDdat$INTENSIVE[match(n_med$MaskID, ACCORDdat$ID)]
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$Visit == time[i] & bp$MaskID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == T)])[1]
  sd.intensive[i] <- sd(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == T)])
  ave.standard[i] <- ave(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == F)])[1]
  sd.standard[i] <- sd(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == F)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == T))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == F))
}
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(110, 150), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average SBP (mmHg)", main = "ACCORD external validation", 
     xaxt = "n")
axis(1, at = 0:10 * 12, labels = 0:10 * 12)
par(new = T)
plot(time, ave.standard,
     xlim = range(time), ylim = c(110, 150), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:7) * 12){
  to.select <- which(n_med$Visit == i & n_med$MaskID %in% ID.high.risk)
  med.intensive[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == T)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == F)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -12, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -24, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -24, cex = 1)
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$Visit == time[i] & bp$MaskID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == T)])[1]
  sd.intensive[i] <- sd(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == T)])
  ave.standard[i] <- ave(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == F)])[1]
  sd.standard[i] <- sd(bp$sbp[to.select][which(bp$INTENSIVE[to.select] == F)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == T))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == F))
}
par(new = T)
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(110, 150), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(time, ave.standard,
     xlim = range(time), ylim = c(110, 150), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:7) * 12){
  to.select <- which(n_med$Visit == i & n_med$MaskID %in% ID.low.risk)
  med.intensive[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == T)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == F)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -24, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -24, cex = 1)
#legends
legend(x = 0, y = 150, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 40, y = 150, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

########## plot delta systolic blood pressure agaist time, subgrouped by high/low risk? ##########
#z test function
z.test <- function(mu1, sd1, mu2, sd2){
  p <- pnorm((mu1-mu2)/sqrt(sd1^2 + sd2^2))
  if (p > 0.5){
    return (2 * (1-p))
  }else{
    return (2 * p)
  }
}

#plot ave, standard red, intensive green
pdf(file = "../../writing/SBP_delta_n_med_ACCORD.pdf", 6, 6)
par(mar = c(8, 8, 4, 2))
#load bp
bp <- read.csv("bloodpressure.csv", stringsAsFactors = F)                   #[181991, 5]
bp$Visit[which(bp$Visit == "BLR")] <- 0
bp <- bp[-which(bp$Visit == "EXIT" | is.na(bp$Visit)), ]
bp$SBP_delta <- rep(NA, dim(bp)[1])
ID <- ACCORDdat$ID
#construct delta
for (i in ID){
  to.select <- which(bp$MaskID == i)
  if(length(bp$sbp[to.select][which(bp$Visit[to.select] == 0)]) == 0){
    next
  }
  bp$SBP_delta[to.select] <- bp$sbp[to.select] - bp$sbp[to.select][which(bp$Visit[to.select] == 0)]
}
bp$Visit <- as.numeric(gsub("F", "", bp$Visit))
bp$INTENSIVE <- ACCORDdat$INTENSIVE[match(bp$MaskID, ACCORDdat$ID)]
ID.high.risk <- ACCORDdat[pred, "ID"]
ID.low.risk <- ACCORDdat[-pred, "ID"]
time <- sort(unique(bp$Visit))
time <- c(0:3, (1:21) * 4)
#load n_med
n_med <- read.csv("concomitantmeds.csv", stringsAsFactors = F)               #[60560, 56]
names(n_med)
unique(n_med$Visit)
summary(n_med)
n_med$Visit[which(n_med$Visit == "BLR")] <- 0
n_med$Sum <- rowSums(n_med[, c(3:5, 7:16)], na.rm = T)
n_med <- n_med[-which(n_med$Visit == "EXIT"), ]
n_med$Visit <- as.numeric(gsub("F", "", n_med$Visit))
n_med$INTENSIVE <- ACCORDdat$INTENSIVE[match(n_med$MaskID, ACCORDdat$ID)]
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$Visit == time[i] & bp$MaskID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == T)])[1]
  sd.intensive[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == T)])
  ave.standard[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == F)])[1]
  sd.standard[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == F)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == T))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == F))
}
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(-30, 0), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average SBP (mmHg)", main = "ACCORD external validation", 
     xaxt = "n")
axis(1, at = 0:10 * 12, labels = 0:10 * 12)
par(new = T)
plot(time, ave.standard,
     xlim = range(time), ylim = c(-30, 0), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:7) * 12){
  to.select <- which(n_med$Visit == i & n_med$MaskID %in% ID.high.risk)
  med.intensive[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == T)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == F)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -12, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -24, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -24, cex = 1)
write.csv(file = "../../writing/temp.csv", data.frame(ave.intensive, sd.intensive/sqrt(n.intensive), ave.standard, sd.standard/sqrt(n.standard)))
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$Visit == time[i] & bp$MaskID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == T)])[1]
  sd.intensive[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == T)])
  ave.standard[i] <- mean(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == F)], na.rm = T)[1]
  sd.standard[i] <- sd(bp$SBP_delta[to.select][which(bp$INTENSIVE[to.select] == F)], na.rm = T)
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == T))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == F))
}
write.csv(file = "../../writing/temp.csv", data.frame(ave.intensive, sd.intensive/sqrt(n.intensive), ave.standard, sd.standard/sqrt(n.standard)))
par(new = T)
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(-30, 0), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(time, ave.standard,
     xlim = range(time), ylim = c(-30, 0), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:7) * 12){
  to.select <- which(n_med$Visit == i & n_med$MaskID %in% ID.low.risk)
  med.intensive[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == T)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == F)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -24, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -24, cex = 1)
#legends
legend(x = 0, y = -25, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 40, y = -25, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

########## plot diastolic blood pressure agaist time, subgrouped by high/low risk? ##########
#plot ave, standard red, intensive green
pdf(file = "../../writing/DBP_n_med_ACCORD.pdf", 6, 6)
par(mar = c(8, 8, 4, 2))
#load bp
bp <- read.csv("bloodpressure.csv", stringsAsFactors = F)                   #[181991, 5]
bp$Visit[which(bp$Visit == "BLR")] <- 0
bp <- bp[-which(bp$Visit == "EXIT" | is.na(bp$Visit)), ]
bp$Visit <- as.numeric(gsub("F", "", bp$Visit))
bp$INTENSIVE <- ACCORDdat$INTENSIVE[match(bp$MaskID, ACCORDdat$ID)]
ID.high.risk <- ACCORDdat[pred, "ID"]
ID.low.risk <- ACCORDdat[-pred, "ID"]
time <- sort(unique(bp$Visit))
time <- c(0:3, (1:21) * 4)
#load n_med
n_med <- read.csv("concomitantmeds.csv", stringsAsFactors = F)               #[60560, 56]
names(n_med)
unique(n_med$Visit)
summary(n_med)
n_med$Visit[which(n_med$Visit == "BLR")] <- 0
n_med$Sum <- rowSums(n_med[, 3:56], na.rm = T)
n_med <- n_med[-which(n_med$Visit == "EXIT"), ]
n_med$Visit <- as.numeric(gsub("F", "", n_med$Visit))
n_med$INTENSIVE <- ACCORDdat$INTENSIVE[match(n_med$MaskID, ACCORDdat$ID)]
###among high risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$Visit == time[i] & bp$MaskID %in% ID.high.risk)
  ave.intensive[i] <- ave(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == T)])[1]
  sd.intensive[i] <- sd(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == T)])
  ave.standard[i] <- ave(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == F)])[1]
  sd.standard[i] <- sd(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == F)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == T))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == F))
}
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(60, 100), type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time (Month)", ylab = "Average DBP (mmHg)", main = "ACCORD external validation", 
     xaxt = "n")
axis(1, at = 0:10 * 12, labels = 0:10 * 12)
par(new = T)
plot(time, ave.standard,
     xlim = range(time), ylim = c(60, 100), type = "l", col = "red", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:7) * 12){
  to.select <- which(n_med$Visit == i & n_med$MaskID %in% ID.high.risk)
  med.intensive[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == T)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 3.8, at = i, cex = 1)
  med.standard[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == F)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 4.6, at = i, cex = 1)
  j <- j + 1
}
mtext("Mean number of medications", side = 1, line = 3.2, at = -12, cex = 1)
mtext("high risk, intensive", side = 1, line = 3.8, at = -24, cex = 1)
mtext("high risk, standard", side = 1, line = 4.6, at = -24, cex = 1)
###among low risk patients
ave.intensive <- c()
ave.standard <- c()
sd.intensive <- c()
sd.standard <- c()
n.intensive <- c()
n.standard <- c()
for(i in 1:length(time)){
  to.select <- which(bp$Visit == time[i] & bp$MaskID %in% ID.low.risk)
  ave.intensive[i] <- ave(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == T)])[1]
  sd.intensive[i] <- sd(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == T)])
  ave.standard[i] <- ave(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == F)])[1]
  sd.standard[i] <- sd(bp$dbp[to.select][which(bp$INTENSIVE[to.select] == F)])
  n.intensive[i] <- length(which(bp$INTENSIVE[to.select] == T))
  n.standard[i] <- length(which(bp$INTENSIVE[to.select] == F))
}
par(new = T)
plot(time, ave.intensive, 
     xlim = range(time), ylim = c(60, 100), type = "l", col = "green", lwd = 2, lty = 2,
     xlab = "", ylab = "", main = "", xaxt = "n")
par(new = T)
plot(time, ave.standard,
     xlim = range(time), ylim = c(60, 100), type = "l", col = "green", lwd = 2, lty = 1,
     xlab = "", ylab = "", main = "", xaxt = "n")
arrows(time, ave.intensive + sd.intensive/sqrt(n.intensive), time, ave.intensive - sd.intensive/sqrt(n.intensive),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
arrows(time, ave.standard + sd.standard/sqrt(n.standard), time, ave.standard - sd.standard/sqrt(n.standard),
       length = 0.02, angle = 90, code = 3, col = "black", lwd = 2)
#add number of medications
med.intensive <- c()
med.standard <- c()
j <- 1
for (i in (0:7) * 12){
  to.select <- which(n_med$Visit == i & n_med$MaskID %in% ID.low.risk)
  med.intensive[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == T)], na.rm = T)[1], 1)
  mtext(med.intensive[j], side = 1, line = 5.4, at = i, cex = 1)
  med.standard[j] <- round(mean(n_med$Sum[to.select][which(n_med$INTENSIVE[to.select] == F)], na.rm = T)[1], 1)
  mtext(med.standard[j], side = 1, line = 6.2, at = i, cex = 1)
  j <- j + 1
}
mtext("low risk, intensive", side = 1, line = 5.4, at = -24, cex = 1)
mtext("low risk, standard", side = 1, line = 6.2, at = -24, cex = 1)
#legends
legend(x = 0, y = 100, legend = c("low risk, intensive", "low risk, standard"), 
       col = c("green", "green"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
legend(x = 40, y = 100, legend = c("high risk, intensive", "high risk, standard"), 
       col = c("red", "red"), lty = c(2, 1), bty = "n", lwd = 1.5, cex = 0.9)
dev.off()

rm(list = ls())
sessionInfo()