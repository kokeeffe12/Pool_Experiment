require(tidyr)
require(scales)
require(ggplot2)
require(dplyr)
require(stringr)
require(nlme)
require(grid)
require(broom)
require(agricolae)
library(devtools)
require(DiagrammeR)

# ggplot theme----
theme_figs <- theme_bw() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) 
theme_set(theme_figs)

update_geom_defaults('point', c(size = 1.5))
update_geom_defaults('errorbar', c(size = 0.5))

# set up data----
#Prevalence Data
setwd("~/Documents/Mitchell Lab/Brandon/Summer 2018/Pool Experiment/Data")

poolprev.dat <- read.csv("Kayleigh_PoolData.csv")  %>%
  filter(Dam == "D" | Dam == "ND") %>%
  filter(Distance == "12" | Distance == "24" | Distance == "36")

poolprev.dat
summary(poolprev.dat)
poolprev.dat$Pool<-as.factor(poolprev.dat$Pool)
poolprev.dat$fDistance<-as.factor(poolprev.dat$Distance)
poolprev.dat$Distance<-as.numeric(poolprev.dat$Distance)

# restructure data
poolprev.long <- gather(poolprev.dat, key = Parasite, value = Infected, Rhiz, Coll, Pucc, Sept) %>% 
  mutate(Parasite = factor(case_when(.$Parasite == "Coll" ~ "Colletotrichum",
                                     .$Parasite == "Rhiz"~ "Rhizoctonia",
                                     .$Parasite == "Pucc" ~ "Puccinia",
                                     .$Parasite == "Sept" ~ "Septoria"), levels = c("Colletotrichum", "Rhizoctonia", "Puccinia", "Septoria")))
# plot prevalence survey data
ss.pp <- ggplot(poolprev.long[poolprev.long$Parasite=="Rhizoctonia",], aes(x = DAE, y = Infected, group = Pool, color=Pool_Endo)) + 
  facet_wrap(~Parasite, scales = "free")+ 
  stat_summary(geom = "point", fun.y = "mean", size = 10, shape = "_") + 
  stat_summary(geom = "line", fun.y = "mean", size = 1) + 
  xlab(label = "Proportion of Leaves Infected with Rhizoctonia")+
  theme(legend.position = c(.3,.85)) 
ss.pp

#Restructure again to summarize the pool data 
poolprev.long2<-poolprev.long[poolprev.long$Parasite=="Rhizoctonia",]

#Pool_Endo only
PP.sum<-poolprev.long2 %>%
  group_by(Pool, DAE, Pool_Endo) %>%
  summarize(Rhiz.Prev = mean(Infected, na.rm = T)) %>%
  ungroup()

Rmod.null <-gls(Rhiz.Prev~1, method="REML", data=PP.sum)
Rmod.gls <-gls(Rhiz.Prev~Pool_Endo*DAE, method="REML", data=PP.sum)

Rmod.ri<-lme(Rhiz.Prev~Pool_Endo*DAE, 
             random=~1|Pool, method="REML", 
             data=PP.sum)

Rmod.slopes<-lme(Rhiz.Prev~Pool_Endo*DAE, 
                 random=~DAE|Pool, method="REML", 
                 control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
                 data=PP.sum)

AIC(Rmod.null, Rmod.gls, Rmod.ri, Rmod.slopes)
anova(Rmod.gls, Rmod.ri, Rmod.slopes)
#random intercepts has lowest AIC, but not significant with random slopes
#random slopes takes into account repeated measures so keeping
resids.fig(Rmod.slopes, PP.sum) #some heteroscadasticity?

#Add in temporal autocorrelation structure
Rmod.slopes2<-lme(Rhiz.Prev~Pool_Endo*DAE, 
             random=~DAE|Pool, method="REML", 
             correlation=corCAR1(form=~DAE|Pool),
             data=PP.sum)

#Add in third-order polynomial
Rmod.slopes3<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 3, raw=TRUE), 
              random=~DAE|Pool, method="REML", 
              correlation=corCAR1(form=~DAE|Pool),
              data=PP.sum)

Rmod.slopes4<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 3, raw=TRUE), 
              random=~DAE|Pool, method="REML", 
              data=PP.sum)

AIC(Rmod.slopes, Rmod.slopes2)
anova(Rmod.slopes, Rmod.slopes2) #correlation structure is significantly better
AIC(Rmod.slopes3, Rmod.slopes4)
anova(Rmod.slopes3, Rmod.slopes4) #correlation structure is significantly better

resids.fig(Rmod.slopes3, PP.sum) #better, homoscedastic

Rmod.ml<-update(Rmod.slopes3, method="ML")
require(car)
anova(Rmod.ml)
#no interaction
Rmod2 <- update(Rmod.ml, .~. -Pool_Endo:poly(DAE, 3, raw = TRUE))
anova(Rmod.ml, Rmod2)
anova(Rmod2)

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Rmod2 %>% broom::augment() %>% 
  ggplot(aes(x=DAE,y=Rhiz.Prev))+
  theme_bw() + 
  #geom_line(data =  PP.sum, aes(x=DAE, y=Rhiz.Prev, group = Plant, color=Pool_Endo, alpha=0.00001))+
  #geom_smooth(aes(y=.fixed),size=2, method = "lm", se = F)+ #this is what you want to change.
  geom_point(aes(color=Pool_Endo),size=3, position=position_jitter(width=0.5), alpha=0.3)+
  geom_line(aes(y=.fixed, color=Pool_Endo),size=3)+
  scale_color_manual(values=cbPalette)+
  ylab("Proportion of Leaves \nExhibiting Rhizoctonia Lesions") +
  xlab("Days after Rhizoctonia Inoculation")+
  guides(color=guide_legend(title="Endophyte Inoculation"))+
  theme(axis.text=element_text(size=12),
        axis.title= element_text(size=12),
        axis.ticks=element_line(size=2),
        legend.title = element_text(size=12),
        legend.position=c(0.27, 0.83),
        legend.text=element_text(size=12),
        legend.background = element_blank())+
  scale_color_manual(name="Endophyte Inoculation", 
                          labels = c("Absent", 
                                     "Present"), 
                          values = c("Endo_Free"="#E69F00", 
                                     "Endo_Inoc"="#56B4E9"))

#Maximum prevalence
MaxPP<-PP.sum %>%
  group_by(Pool, Pool_Endo) %>%
  summarize(MaxRhiz.Prev = max(Rhiz.Prev, na.rm = T)) %>%
  ungroup()

t.test(MaxRhiz.Prev~Pool_Endo, data=MaxPP)

require(plotrix)
myData <- aggregate(MaxPP$MaxRhiz.Prev,
                    by = list(cyl = MaxPP$Pool_Endo),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n = length(x)))
myData <- do.call(data.frame, myData)
myData$se <- myData$x.sd / sqrt(myData$x.n)

colnames(myData) <- c("Pool_Endo", "mean", "sd", "n", "se")

myData$names <- c(paste(myData$Pool_Endo, "Endophyte Inoculation"))

ggplot(myData)+
  geom_bar(aes(x=Pool_Endo, y=mean, fill=Pool_Endo), width=0.6, stat="identity", alpha=0.7)+
  geom_errorbar(aes(x=Pool_Endo, ymin=mean-se, ymax=mean+se), alpha=0.9, size=1.1, width=0.3)+
  theme(legend.position="none")+
  scale_fill_manual(name="Endophyte Inoculation", 
                     labels = c("Absent", 
                                "Present"), 
                     values = c("Endo_Free"="#E69F00", 
                                "Endo_Inoc"="#56B4E9"))+
  ylab("Peak Rhizoctonia Prevalence") +
  xlab("Endophyte Inoculation ")+
  scale_x_discrete(breaks=c("Endo_Free","Endo_Inoc"),
                   labels=c("Absent", "Present"))

resids.fig <- function(mod, df) {
  residdf <- dplyr::mutate(df, resids = residuals(mod, type = 'normalized'),
                           fits = fitted(mod))
  fig2 <-ggplot(residdf, aes(x = fits, y = resids)) + geom_point() +
    labs(x = 'Fitted values', y = '')
  
  fig3 <- ggplot(residdf) + stat_qq(aes(sample = resids)) +
    labs(x = 'Theoretical Quantiles', y = 'Sample Quantiles')
  
  # qqline plot = FALSE, according to James should work
  
  fig4 <- ggplot(residdf, aes(x = resids)) + geom_histogram(aes(y=..density..), colour = 'grey50') +
    labs(x = 'Residuals', y = 'Frequency') + scale_y_continuous(expand = c(0, 0)) +
    stat_function(fun = dnorm, color = "red", args = list(mean = mean(residdf$resids),
                                                          sd = sd(residdf$resids)))
  grid.draw(rbind(ggplotGrob(fig2), ggplotGrob(fig3), ggplotGrob(fig4), size = 'first'))
  
  return(summary(mod))
}
