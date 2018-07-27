### R code from vignette source '/home/james/work/stats_consult/pkelly/powernew.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: powernew.Rnw:12-14
###################################################
library(lme4)
library(ggplot2)


###################################################
### code chunk number 2: dd
###################################################
opts_chunk$set(warning=FALSE) 


###################################################
### code chunk number 3: powernew.Rnw:33-53
###################################################
#n is the number of leaves surveyed per plant
#baserhiz=rhizoctonia prevalence at center (not sure if this should be of time or endophyte)
#timenum=number of days
#poolnum=total number of pools
#endonum=total number of endo treatments
#plant num=number of plants surround central
#NOTE: still haven't included distance
rhizsimplot <-function(n=20, 
                       baserhiz=0.2, 
                       endoeffect=-.1, 
                       timeeffect=0.02, 
                       timenum=21,
                       poolnum=28,
                       endonum=2, 
                       plantnum=12,
                       poolvar=.001, 
                       endovar=.02){
  
  require(ggplot2)
  endo=c(rep("0", poolnum/2), rep("1", poolnum/2)) #half have endo, half don't
  plant <- factor(1:plantnum)
  pool <- factor(1:poolnum)
  time<-as.numeric(1:timenum)
  
  
  data<-data.frame(
    endo=as.numeric(endo),
    pool=rep(pool, each=plantnum),
    plant=rep(plant, poolnum)
  )
  data2<-do.call("rbind", replicate(timenum, data, simplify = FALSE))#replicate original dataframe for each day of surveying
  data2$time<-rep(time, each=nrow(data))
  data2$centendo<-data2$endo-mean(data2$endo) #not sure if this fits here since this is a factor 
  data2$centtime<-as.numeric(data2$time)-mean(data2$time) #following the previous example
  pooleffect <- rnorm(poolnum,0,poolvar)
  nu <- with(data2,log(baserhiz/(1-baserhiz))+endo*endoeffect+time*timeeffect+pooleffect[as.numeric(pool)]) #stumped if this is correct
  data2$rhiz <- sapply(nu,function(x) rbinom(1,size=n,prob=exp(x)/(1+exp(x))))
  data2$healthy <- n-data2$rhiz
  ggplot(data2,aes(y=rhiz/(rhiz+healthy),x=time,color=endo))+geom_point()+geom_line(aes(y=exp(log(baserhiz/(1-baserhiz))+endoeffect*endo)/(1+exp(log(baserhiz/(1-baserhiz))+endoeffect*centendo))))
}

rhizsimplot()


###################################################
### code chunk number 4: powernew.Rnw:58-79
###################################################

rhizsimfit <- function(n=20, baserhiz=0.2, endoeffect=.05, disteffect=0.05, distnum=3, poolnum=28, endonum=2, plantnum=12, poolvar=.0001, endovar=.02){
  endo=c(rep("0", poolnum/2), rep("1", poolnum/2))
  plant <- factor(1:plantnum)
  pool <- factor(1:poolnum)
  dist<-as.numeric(rep(1:distnum, each=4))
  
  data<-data.frame(
    endo=as.numeric(endo),
    pool=rep(pool, each=plantnum),
    dist=rep(dist, poolnum),
    plant=rep(plant, poolnum)
  )
  pooleffect <- rnorm(poolnum,0,poolvar)
  data$centdist<-data$dist-mean(data$dist)
  data$centendo<-data$endo-mean(data$endo)
  nu <- with(data,log(baserhiz/(1-baserhiz))+endo*endoeffect+dist*disteffect+pooleffect[as.numeric(pool)])
  data$rhiz <- sapply(nu,function(x) rbinom(1,size=n,prob=exp(x)/(1+exp(x))))
  data$healthy <- n-data$rhiz
  trymod <- glmer(cbind(rhiz,healthy)~endo+dist+(1|pool),data=data,family=binomial)
  trymod2 <- glmer(cbind(rhiz,healthy)~+(1|pool),data=data,family=binomial)
  x <- anova(trymod,trymod2)$"Pr(>Chisq)"[2]
  x
}

sum(replicate(100, rhizsimfit()<0.05))/100
> 