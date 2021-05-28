#Load Packages
library(metafor)
library(dplyr)

#Import Data & calculate effect and variance
htn <- read.csv("HTN.csv",  header = TRUE)
htn <- htn[1:11,] %>% #subset cross-sectional studies only
  mutate(yi=log(Odds.ratio),
         vi=(log(Upper.Limit)-log(Lower.Limit))/(2*1.96)
         ^2)

#Run random effects model
htnRes <- rma(yi, vi, data=htn, method='REML')
htnRes
confint(htnRes)

#Back-transform results for easier interpretation
htnResOR <- predict(htnRes, transf=exp, digits=3)
htnResOR

#Create Forest plot
tiff(filename="HTN Forest Plot Delta Hat.tif", width=16.80104, height=12.35604, units='cm', res=300)
forest(htn$Odds.ratio, #use odds ratios for interpretability
       ci.lb=htn$Lower.Limit,
       ci.ub=htn$Upper.Limit,
       rows=13:3, #to leave space for summary effect and abline
       digits=c(2,1),
       cex=0.8,
       xlab="Odds Ratio",
       xlim=c(-1,2.8),
       alim=c(0,2),
       ylim=c(0,15.8), #to get top abline in correct position
       slab=htn$Study.name,
       refline=1 #reference line to 1 instead of 0
)
addpoly(htnRes, row=1, transf=exp, cex=0.8) #summary effect
abline(h=2) #horizontal line below final study

#add annotations       
text(-1, 14.3, "Authors, Year", pos=4, cex=0.8)
text(2.8, 14.3, "OR [95% CI]", pos=2, cex=0.8)

dev.off()