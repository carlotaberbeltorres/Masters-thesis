## setting the working directory.

setwd("C:/Users/cb821/Desktop/project")

## Loading the size scan files into my environment

scans<- read.csv("scans.csv")


## Loading the mortality data 

treatments <- read.csv("microbes.csv")

## downloading packages I might need for survival analysis

#install.packages("survival")
library(survival)
#install.packages("survminer")
library(survminer)
#install.packages("lubridate")
library(lubridate)
#install.packages("ranger")
library(ranger)
library(ggplot2)
library(dplyr)
#install.packages("ggfortify")
library(ggfortify)

## deleting the water and feed treatment
experiment <- filter(treatments, Treatments !="Control")
experiment

## merge datasets

Final <- merge(scans, treatments, by="LarvaeID")

## Looking at survival curves per treatment only taking B+ out to see the effect of pesticides on larval survival

# Let's exclude B+

Pesticides.BMinus <- filter(treatments, treatment == "A-Untreated")

## deleting the water and feed treatment

Pest.BMinus <- filter(Pesticides.BMinus, Treatments !="Control")
##Survival analysis
Pesticides.Survival <- survfit(Surv(DOD, cencorDOD) ~ Treatments, data=Pest.BMinus)

summary(Pesticides.Survival)

## Plot survival with pesticide exposure data
ex<-ggsurvplot(
  fit = survfit(Surv(DOD, cencorDOD) ~ Treatments, data = Pest.BMinus), 
  xlab = "Time (Days)", 
  ylab = "Overall survival probability",
  legend.labs = c("Control","Thiacloprid","Lamda-Cyhalothrin","Glyphosate","Azoxystrobin"), 
  legend=c("bottom") ,font.legend=12, legend.title = "",
  pval = TRUE, pval.size= 5,pval.coord = c(0.1, 0.1),
  title = "Larval survival under pesticide exposure",font.main = c(16, "bold"))
ex




## Calculating cox regression

Cox.Pesticides <- coxph(Surv(DOD, cencorDOD) ~Treatments  +colony_id,  data = Pest.BMinus)
summary(Cox.Pesticides)
tab_model(Cox.Pesticides)

## Forest plot
ggforest.Pesticides <- ggforest(Cox.Pesticides, data = Pest.BMinus, cpositions=c(0.01,0.15,0.3), 
                          main = "Hazard ratios for B+ and B-", fontsize = 1)
ggforest(Cox.Pesticides)

## Looking at whether B+ affects larval survival 

Microbes <- Final %>% filter(Treatment %in% c("Acetone", "1"))

Microbes.Survival <- survfit(Surv(DOD, cencorDOD) ~ Treatment, data=Microbes)
autoplot(Microbes.Survival)

## Calculating cox regression

Cox.Microbes <- coxph(Surv(DOD, cencorDOD) ~ Treatment+treatment ,  data = Microbes)
summary(Cox.Microbes)

tab_model(Cox.Microbes)

##Forest Plot

ggforest.Microbes <- ggforest(Cox.Microbes, data = Microbes, cpositions=c(0.01,0.15,0.3), 
                          main = "Hazard ratios for B+ and B-", fontsize = 1)
ggforest(Cox.Microbes)

##Looking at pesticides and microbes

All.Treat <- filter(treatments, Treatments !="Control")
All <- survfit(Surv(DOD, cencorDOD) ~ Treatments, data=All.Treat)
autoplot(All)

## Calculating cox regression

Cox.All <- coxph(Surv(DOD, cencorDOD) ~Treatments + treatment + colony_id ,  data = All.Treat)
summary(Cox.All)

#plot model output into a table
library(sjPlot)

tab_model(Cox.All)
## Forest plots

ggforest.data <- ggforest(Cox.All, data = All.Treat, cpositions=c(0.01,0.15,0.3), 
                                main = "Hazard ratios for B+ and B-", fontsize = 1)
ggforest(Cox.All)


tab_model(Cox.Pesticides, Cox.Microbes,Cox.All, CSS = list(css.depvarhead = '+color: blue;'), file="Cox.doc")
### Lets analyse the change in size
###Calculating larval growth rate
growth <- Final %>%
  arrange(LarvaeID,Treatments,colony_id,treatment) %>%
  group_by(LarvaeID, Treatments,treatment) %>% 
  mutate(base_rate = first(Area.Start),
         growth_rate = ((Area.Finish - base_rate)/DOD*100) )%>%
  select(-base_rate)

head(growth)
scans <- filter(Final, Treatments !="Control")
scans

##load package to do GLMs
library(lmerTest)

## Model looking only at pesticides
M1<- filter(growth, treatment == "A-Untreated")
M1treat<-filter(M1, Treatments !="Control")
m1<- glm((growth_rate)~(Treatments)+(colony_id), data = M1treat)
summary(m1)


##run again to reset variables
growth <- Final %>%
  arrange(LarvaeID,Treatments,colony_id,treatment) %>%
  group_by(LarvaeID, Treatments,treatment) %>% 
  mutate(base_rate = first(Area.Start),
         growth_rate = ((Area.Finish - base_rate)/DOD*100) )%>%
  select(-base_rate)

##Model looking at only B+
M2<- growth %>% filter(Treatment %in% c("Acetone", "1"))
m2 <- glm((growth_rate)~(Treatment)+(treatment)+(colony_id), data=M2)
summary(m2)

##run again to reset variables
growth <- Final %>%
  arrange(LarvaeID,Treatments,colony_id,treatment) %>%
  group_by(LarvaeID, Treatments,treatment) %>% 
  mutate(base_rate = first(Area.Start),
         growth_rate = ((Area.Finish - base_rate)/DOD*100) )%>%
  select(-base_rate)


## Model considering both influence of pesticides and microbes
growth2 <- filter(growth, Treatments !="Control")
m3<- glm((growth_rate)~ (Treatments)+(treatment)+(colony_id) ,data=growth2)
summary(m3)
plot(m3)

## combining model outputs onto a table
library(sjPlot)
tab_model(m1,m2, m3)
tab_model(m1, m2,m3, CSS = list(css.depvarhead = '+color: blue;'), file="m1m2.doc")

library(ggplot2)

growth <- filter(growth, Treatment.x !="1")



growth <- filter(growth, Treatments !="Control")
growth_BMinus <- filter(growth, treatment == "A-Untreated")
BMinus<-ggplot(growth_BMinus, aes(x=Treatments, y=growth_rate, fill=colony_id)) + 
  geom_boxplot()+labs(x = "Treatments", y = "Larval growth (mm)",
                                                title = "Untreated",font.main = c(16, "bold"))

BMinuss<- BMinus + scale_y_continuous(limits=c(0,6000),breaks = seq(0, 6000,1000))


growth_BPlus <- filter(growth, treatment == "Supplemented")
BPlus<-ggplot(growth_BPlus, aes(x=Treatments, y=growth_rate, fill=colony_id)) + 
  geom_boxplot()+labs(x = "Treatments", y = "Larval growth (mm)",
                                                title = "Supplemented",font.main = c(16, "bold"),ylim(0,8000))
BPluss<-BPlus+ scale_y_continuous(limits=c(0,6000),breaks = seq(0, 6000,1000))

## stacking both graphs together
figure <- ggarrange(BPluss, BMinuss,
                 ncol = 1, nrow = 2)
figure



##References needed for my project 
citation()
citation("survival")             
citation ("lmerTest")
