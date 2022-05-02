## LIBRARY
library(stats4,maxLik,rootSolve)
library(tidyverse)

###################
## ASN FUNCTIONS ##
###################
# ASN DENSITY (PDF)
dASN<-function(t,mu,sigma,alpha){
  den<-( ((1-alpha*((t-mu)/sigma))^2 + 1)/(sigma*(2+alpha^2)))*dnorm((t-mu)/sigma)
  return(den)
}

# ASN CUMULATIVE (CDF)
pASN<-function(t,mu,sigma,alpha) {
  cum<- pnorm((t-mu)/sigma) + alpha*( (2*sigma - alpha*(t-mu)) / (sigma*(2+alpha^2)) ) * dnorm((t-mu)/sigma)
  return(cum)
}

# ASN RANDOM NUMBER GENERATOR
rASN<-function(n,pmu,psigma,palpha){
  t<-c() 
  a<-runif(n,0,1)
  lim1<-100
  for (i in 1:N) 
    t[i]<-uniroot(acharxp, lower = -100, upper = 100, mu=pmu,sigma=psigma,alpha=palpha,u=a[i])$root 
  return(t)
}

acharxp<-function(t,mu,sigma,alpha,u){
  u-pASN(t,mu,sigma,alpha)
}

## ASN Log-LIKELIHOOD
loglike<-function(theta) {
  mu<-theta[1]
  sigma<-theta[2]
  alpha<-theta[3]
  den<-dASN(t,mu,sigma,alpha)
  aux<-sum(log((1-alpha*((t-mu)/sigma))^2 + 1))-n*log(sigma)-n*log((2+alpha^2))+sum(dnorm((t-mu)/sigma, log = TRUE))
  return(aux)
} 

## ASN RIGHT-Tail ANDERSON-DARLING (RADE)
frmad<- function(theta){
  mi<-theta[1]
  sigma<-theta[2]
  alfa<-theta[3]
  auxh<- -n/2+2*sum(pASN(t,mi,sigma,alfa))+(1/n)*sum((2*seq(1,n,1)-1)*log(1-pASN(taux,mi,sigma,alfa)))
  return(auxh)
}  

## ASN Maximum Product of Spacings (MPS)
fmps<- function(theta){
  mi<-theta[1]
  sigma<-theta[2]
  alfa<-theta[3]
  t<-sort(dados, decreasing = F)
  n<-length(dados)-1
  D<-NULL
  D[1]<-pASN(t[1],mi,sigma,alfa)
  D[n+1]<-1
  for (i in 2:n)
  {
    if(t[i] == t[i-1])
    {
      D[i] = dASN(t[i],mi,sigma,alfa)
      next
    } 
    
    D[i]<-pASN(t[i],mi,sigma,alfa)-pASN(t[i-1],mi,sigma,alfa) }
  
  aux<-1/(n+1)*sum(log(D))
  return(aux)
} 


#################################
##  DESCRIPTIVE
DataSet=read.csv("https://raw.githubusercontent.com/ProfNascimento/ASN/main/AtacamaRivers_Water.csv",header = T)

DB = DataSet %>%
  gather("X2011","X2012","X2013","X2014","X2015","X2016","X2017","X2018","X2019","X2020", key = Year, value = Flux)

ggplot(DB, aes(x = Flux, y = ..density..)) +
  geom_histogram(alpha = 0.3, bins = 100) +
  ylab("Density")+
  xlab("Monthly Flux")+
  geom_density(size = .5, color = "red")

DB$MONTH = factor(DB$MONTH,levels=c("ENE","FEB","MAR","ABR","MAY","JUN","JUL","AGO","SEP","OCT","NOV","DIC"))
levels(DB$MONTH)=c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

tapply(DB$Flux, DB$MONTH, summary)

ggplot(DB, aes(MONTH, log(Flux))) +
  geom_smooth(aes(y=Flux,x=MONTH), data=DB, method = "loess", span = 0.75)+
  geom_jitter(width = 0.1, height = 0.1) +
  facet_wrap(~Year) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dados=log(as.numeric(unlist(c(na.omit(DataSet[,3:12])))))
hist(dados ,breaks=20 ,freq=FALSE, xlab="Flux", main = "Histogram - Monthly log(Water Flux)")

######################
## ESTIMATION STEPS ##
######################
old=options(digits=4)
set.seed(123456)

t<-sort(dados)
n<-length(dados)
taux<-sort(dados, decreasing = TRUE)

# Initial Mu 
pmu<-mean(t)
# Initial Sigma
psigma<-sd(t)
# Alpha grid search
w<-seq(-10,10,0.5)
aux1<-c()
for(i in 1:length(w)){
  ESThat  <- try(maxBFGS(loglike, start=c(pmu,psigma,w[i]))$estimate)
  aux1[i]<-ks.test(t,pASN,ESThat[1],ESThat[2],ESThat[3])$p.value
}
(palpha<- w[which.max(aux1)])
(ASNmle  <- try(maxBFGS(loglike, start=c(pmu,psigma,palpha))$estimate) )

aux2<-c()
for(i in 1:length(w)){
  ESThat  <- try(maxBFGS(frmad, start=c(pmu,psigma,w[i]))$estimate)
  aux2[i]<-ks.test(t,pASN,ESThat[1],ESThat[2],ESThat[3])$p.value
}
(palpha<- w[which.max(aux2)])
(ASNmadr  <- try(maxBFGS(frmad, start=c(pmu,psigma,palpha))$estimate) )

aux3<-c()
for(i in 1:length(w)){
  ESThat  <- try(maxBFGS(frmad, start=c(pmu,psigma,w[i]))$estimate)
  aux3[i]<-ks.test(t,pASN,ESThat[1],ESThat[2],ESThat[3])$p.value
}
(palpha<- w[which.max(aux3)])
(ASNmps  <- try(maxBFGS(fmps, start=c(pmu,psigma,palpha))$estimate) )


## ASN Kolmogorov-Smirnov test
ks.test(t,pASN,ASNmle[1],ASNmle[2],ASNmle[3])
ks.test(t,pASN,ASNmadr[1],ASNmadr[2],ASNmadr[3])
ks.test(t,pASN,ASNmps[1],ASNmps[2],ASNmps[3])

## ASN Akaike's Information Criterion (AIC)
-2*loglike(c(ASNmle[1],ASNmle[2],ASNmle[3]))+6
-2*loglike(c(ASNmadr[1],ASNmadr[2],ASNmadr[3]))+6
-2*loglike(c(ASNmps[1],ASNmps[2],ASNmps[3]))+6


## VISUAL GOODNESS-OF-FIT
Xhat<-seq(min(t),max(t),0.01)
Ymle<-dASN(Xhat,ASNmle[1], ASNmle[2], ASNmle[3])
Ymadr<-dASN(Xhat,ASNmadr[1], ASNmadr[2], ASNmadr[3])
Ymps<-dASN(Xhat,ASNmps[1], ASNmps[2], ASNmps[3])

par(mfrow=c(1,2))
hist(t,freq =FALSE,ylim = c(0,0.35), main="",xlab="Water Flux (natural log)")
points(Xhat,Ymle,type="l",col="1",lty=2)
points(Xhat,Ymadr,type="l",col="2",lty=2)
points(Xhat,Ymps,type="l",col="3",lty=2)
legend("topright",c("MLE","rADE","MPS"),col=c(1,2,3),lty=c(2,2,2))

#Adjusting KM
ekm<-survival::survfit(Surv(t,rep(1,n))~1)
time<-ekm$time; st<-ekm$surv

ekm<-ecdf(t)
Cmle <- pASN(t,ASNmle[1],ASNmle[2],ASNmle[3])
Cmadr <- pASN(t,ASNmadr[1],ASNmadr[2],ASNmadr[3])
Cmps <- pASN(t,ASNmps[1],ASNmps[2],ASNmps[3])

plot(ekm, xlim=c(min(t),max(t)),ylim=c(0,1),xlab="Monthly log(Water Flux)", ylab="F(X)", do.points = FALSE)
lines(c(min(t),time),c(0,unique(Cmle)), col="1",lty=1)
lines(c(min(t),time),c(0,unique(Cmadr)), col="2",lty=1)
lines(c(min(t),time),c(0,unique(Cmps)), col="3",lty=1)
