library(survival)
library(Hmisc)
library(survivalROC)
library(risksetROC)
library(dplyr)

load("C:/Users/liangcj/Dropbox/Conferences/Karolinska2022 tutorial/Tutorial/actt.RData")

# id: patient id
# trt: placebo or remdesivir arm
# age: age in years
# news: National Early Warning Score (NEWS)
# time.r: time to recovery
# status.r: 0 = censored; 1 = recovery observed
# time.m: time to mortality
# status.m: 0 = censored; 1 = mortality observed
# time.v: time to ventilation
# status.v: 0 = censored; 1 = ventilation observed
# NOTE: for recovery endpoint, if patient died, they were assumed to have
#       a recovery time of "infinity" and thus were censored at day 28 (study
#       follow-up was 28 days)


# split into placebo and remdesivir datasets
acttp <- actt %>% filter(trt=="Placebo")
acttr <- actt %>% filter(trt=="Remdesivir")


# Descriptive plots
set.seed(65)
jitter1 <- rnorm(516, mean=0, sd=0.1)
jitter2 <- rnorm(516, mean=0, sd=0.1)

with(acttp, plot(time.r + jitter1, news + jitter2, 
                 pch=ifelse(status.r==0, 1, 20), 
                 col=ifelse(status.r==0, rgb(1,0,0,0.25), rgb(0,0,0,0.25)), 
                 main="Placebo (n=512)",
                 xlab="days to recovery", ylab="NEWS", cex=1, xlim=c(0,27.5)))

jitter3 <- rnorm(532, mean=0, sd=0.1)
jitter4 <- rnorm(532, mean=0, sd=0.1)
with(acttr, plot(time.r + jitter3, news + jitter4,
                 pch=ifelse(status.r==0, 1, 20), 
                 col=ifelse(status.r==0, rgb(1,0,0,0.25), rgb(0,0,0,0.25)), 
                 main="Remdesivir (n=531)",
                 xlab="days to recovery", ylab="NEWS", cex=1, xlim=c(0,27.5)))


# Exercise 1: manually code c-index using placebo arm, NEWS, days to recovery
#             compare to results from Hmisc::rcorr.cens()
numerator <- 0
denominator <- 0

for(i in 1:515){
  if(is.na(acttp$news[i])) next # skip missing NEWS observations
  for(j in i:516){
    if(is.na(acttp$news[j])) next # skip missing NEWS observations
    if((acttp$time.r[i] == acttp$time.r[j]) & (acttp$status.r[i] != acttp$status.r[j])){
      # tied times but one event and one censor
      denominator <- denominator + 1
      if(acttp$status.r[i]==1 & (acttp$news[i]<acttp$news[j])){
        numerator <- numerator + 1
      }
      if(acttp$status.r[j]==1 & (acttp$news[j]<acttp$news[i])){
        numerator <- numerator + 1
      }
      if(acttp$news[i] == acttp$news[j]){
        # for ties, add 0.5
        numerator <- numerator + 0.5
      }
    }
    if(acttp$time.r[i] != acttp$time.r[j]){
      if( ((acttp$time.r[i] < acttp$time.r[j]) & acttp$status.r[i]==1) |
          ((acttp$time.r[i] > acttp$time.r[j]) & acttp$status.r[j]==1) 
      ){
        # times not tied and smaller time is an event
        denominator <- denominator + 1
        if((acttp$status.r[i]==1) & (acttp$time.r[i] < acttp$time.r[j]) & (acttp$news[i]<acttp$news[j])){
          # observation i has event, and happens before j's, and has lower NEWS than other observation
          numerator <- numerator + 1
        }
        if((acttp$status.r[j]==1) & (acttp$time.r[j] < acttp$time.r[i]) & (acttp$news[j]<acttp$news[i])){
          # observation j has event, and happens before i's, and has lower NEWS than other observation
          numerator <- numerator + 1
        }
        if(acttp$news[i] == acttp$news[j]){
          # for ties, add 0.5
          numerator <- numerator + 0.5
        }
      }
    }
  }
}
numerator/denominator

with(acttp, rcorr.cens(news, Surv(time.r, status.r)))


# Exercise 2: 95% confidence intervals for AUC^C/D(14) using quantile bootstrap 
#             and survivalROC() using placebo arm, recovery endpoint, NEWS
#             NOTE: may need to use negative NEWS instead of NEWS
#                   otherwise will return 1-AUC
ex2 <- with(acttp, survivalROC(Stime=time.r, status=status.r,
            marker=-news, predict.time=14, method="KM"))

ex2$AUC

with(ex2, plot(FP, TP))
with(ex2, lines(FP,TP))

do.one.sroc.boot <- function(m,t,c, time=14, method="KM"){
  # do one bootstrap given marker m, time t, censoring indicator c
  d <- cbind(m, t, c)
  bi <- sample(1:nrow(d), replace=TRUE)
  survivalROC(d[bi,2], d[bi,3], d[bi,1], predict.time=time, method=method)$AUC
}

bootrep <- replicate(1000, with(acttp, do.one.sroc.boot(-news, time.r, status.r, time=14, method="KM")))
quantile(bootrep, c(0.025, 0.975))


# Exercise 3: manually code the weighted mean rank estimator for AUC^I/D(t)
#             use placebo arm, recovery endpoint, NEWS
#             result should be plot
wmrp <- rep(NA, nrow(acttp))
for(i in 1:nrow(acttp)){
  rs <- which(acttp$time.r >= acttp$time.r[i])
  wmrp[i] <- sum(-acttp$news[i] > -acttp$news[rs], na.rm=TRUE)/length(rs)
}

plot(acttp$time.r[acttp$status.r==1], wmrp[acttp$status.r==1], col="blue",
     xlab="days to recovery", ylab="percentile rank", cex=0.75, pch=20, xlim=c(0,27.5))
abline(h=0.5, col="gray")
lines(loess.smooth(acttp$time.r[acttp$status.r==1], wmrp[acttp$status.r==1], span = 0.2), lwd=2)
