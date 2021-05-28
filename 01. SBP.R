#Load Packages
library(metafor)
library(dplyr)

#Import Data & calculate effect and variance
sbp1 <- read.csv("sbp Ind Mean Diff p-value.csv",  header = TRUE)
sbp1 <- sbp1 %>% 
  mutate(yi=diff.in.means,
         vi=(yi/
           (qt(p=p.value/2, #divide p-value by 2 for two-tailed
               df=(con.sample.size+med.sample.size)-2, #degrees of freedom
               lower.tail=F, #to avoid p-values being reversed
               log.p=F))) #do not require logged p-values
           ^2) #square SEM to calculate variance


sbp2 <- read.csv("sbp Ind Mean SD 1TP.csv",  header = TRUE)
sbp2 <- escalc(measure="MD", 
               vtype="HO", #does not assume homoscedasticity
               m2i=con.mean, m1i=med.mean, 
               sd2i=con.sd, sd1i=med.sd, 
               n2i=con.sample.size, n1i=med.sample.size, 
               data=sbp2)


sbp3 <- read.csv("sbp Ind Mean SD.csv",  header = TRUE)
sbp3 <- sbp3 %>% 
  mutate(meanChangeCon = Control.Post.Mean-Control.Pre.Mean,
         meanChangeMed = Med.Diet.Post.Mean-Med.Diet.Pre.Mean,
         sdChangeCon = sqrt((Control.Pre.SD^2 + Control.Post.SD^2) -
                               ((2*Pre.Post.correlation) * 
                               (Control.Pre.SD * Control.Post.SD))),
         sdChangeMed = sqrt((Med.Diet.Pre.SD^2 + Med.Diet.Post.SD^2) -
                               ((2*Pre.Post.correlation) * 
                               (Med.Diet.Pre.SD * Med.Diet.Post.SD)))
         )

sbp3 <- escalc(measure="MD", 
               vtype="HO", #does not assume homoscedasticity
               m2i=meanChangeCon, m1i=meanChangeMed, 
               sd2i=sdChangeCon, sd1i=sdChangeMed, 
               n2i=Control.Sample.size, n1i=Med.Diet.Sample.size, 
               data=sbp3)


sbp4 <- read.csv("sbp Paired MD p-value.csv",  header = TRUE)

sbp4 <- sbp4 %>% 
  mutate(yi=diff.in.means,
         vi=(yi/
               (qt(p=p.value/2, #divide p-value by 2 for two-tailed
                   df=sample.size-1, #degrees of freedom
                   lower.tail=F, #to avoid p-values being reversed
                   log.p=F))) #do not require logged p-values
         ^2) #square SEM to calculate variance


sbp5 <- read.csv("sbp Paired p-value.csv",  header = TRUE)
sbp5 <- sbp5 %>% 
  mutate(yi = post.mean-pre.mean,
         vi=abs((yi/
               (qt(p=p.value/2, #divide p-value by 2 for two-tailed
                   df=sample.size-1, #degrees of freedom
                   lower.tail=F, #to avoid p-values being reversed
                   log.p=F)))) #do not require logged p-values
         ^2) #square SEM to calculate variance


sbp6 <- read.csv("sbp Paired SDdiff.csv",  header = TRUE)
sbp6 <- sbp6 %>% 
  mutate(yi=mean.diff,
         vi=(sd.diff/sqrt(sample.size))
         ^2)


#Retain relevant columns and join to form final dataframe for analysis
sbpFiles <- lapply(list(sbp1,sbp2,sbp3,sbp4,sbp5,sbp6), 
                   subset, 
                   select = c('Study.name', 'Health.status', 'Study.design', 'yi', 'vi'))

data <- as.data.frame(bind_rows(sbpFiles))

mods <- read.csv("Moderators.csv",  header = TRUE)

final <- merge(data, mods, by = "Study.name")
rm(list=setdiff(ls(), "final"))


#Run random effects model
sbpRes <- rma(yi, vi, data=final, method='REML')
sbpRes 
confint(sbpRes)


#Create Forest plot
studyName <- as.character(final$Study.name)
healthStatus <- as.character(final$Health.status)
studyDesign <- as.character(final$Study.design)

tiff(filename="SBP Forest Plot.tif", width=20.69, height=12, units='cm', res=300)
forest(sbpRes,  
       digits=c(2,1), #digits are y axis and x-axis dps
       cex=.8, #cex is font size
       xlab="Mean Difference (mmHg)", #x-axis title
       slab=NA, #remove study labels and add with ilab below for consistency
       xlim=c(-71,35), #horizontal limits of plot region
       alim=c(-20,20), #axis limits
       steps=6, #6 tick marks
       ilab=cbind(studyName, healthStatus, studyDesign), #study info
       ilab.xpos=c(-72,-47.6,-32), #position of study info
       ilab.pos=4 #left aligned
       ) 

  #add annotations       
text(-72, 15.7, "Authors, Year", pos=4, cex=.8)
text(-47.6, 15.7, "Health status", pos=4, cex=.8)
text(-32.5, 15.7, "Study design", pos=4, cex=.8)
text(35, 15.7, "MD [95% CI]", pos=2, cex=.8)

dev.off()



#DIAGNOSTICS
#Leave one out sensitivity analysis
leave1out(sbpRes)

#Baujat plot - Studies that fall to the top right quadrant of the Baujat plot contribute most to overall heterogeneity and the overall result.
tiff(filename="SBP Baujat.tif", width=20.69, height=12, units='cm', res=300)
baujat(sbpRes)
dev.off()

#Diagnostics for outliers and influential cases - asterisk identifies influential cases
inf <- influence(sbpRes)
print(inf)

tiff(filename="SBP Influence.tif", width=20.69, height=12, units='cm', res=300)
plot(inf)
dev.off()

#Funnel plot
tiff(filename="SBP Funnel Plot.tif", width=14, height=12, units='cm', res=300)
funnel(sbpRes, 
       xlab="Mean Difference (mmHg)", 
       digits=c(0,0), 
       hlines=NA)
dev.off()

#Tests for funnel plot
regtest(sbpRes) #p = 0.1177


#SUBGROUP ANALYSIS
subgroupHS <- rma.mv(yi, vi, 
                     mods = ~Health.status, 
                     random = ~ Health.status | Study.name, #random effects MA for each health status subgroup 
                     struct="DIAG", #allows unequal variance between health status subgroups
                     data=final, 
                     digits=3) #number of decimal places of printed results
print(subgroupHS)


subgroupSD <- rma.mv(yi, vi, 
                     mods = ~Study.design, 
                     random = ~ Study.design | Study.name, #random effects MA for each study design subgroup 
                     struct="DIAG", #allows different heteroscedasticity within each study design subgroup - residual heterogeneity variances are heteroscedastic
                     data=final, 
                     digits=3) #number of decimal places of printed results
print(subgroupSD)


#Meta regression - linear
regDuration <- rma(yi, vi, mods = ~Duration, data=final)
print(regDuration) #p=0.2866

regAge <- rma(yi, vi, mods = ~Age, data=final)
print(regAge) #p=0.0799

regAge2 <- rma(yi, vi, mods = ~I(Age^2), data=final)
print(regAge2) #p=0.1126

regAge3 <- rma(yi, vi, mods = ~Age+I(Age^2), data=final)
print(regAge3) #p=0.2127

regBmi <- rma(yi, vi, mods = ~BMI, data=final)
print(regBmi) #p=0.3244

regBp <- rma(yi, vi, mods = ~Baseline.BP, data=final)
print(regBp) #p=0.0005 - remains sig when controlling for other varibles independently or collectively. R2=74.87%

regBp2 <- rma(yi, vi, mods = ~I(Baseline.BP^2), data=final)
print(regBp2) #p=0.0005. R2=74.19%

regBp3 <- rma(yi, vi, mods = ~Baseline.BP+I(Baseline.BP^2), data=final)
print(regBp3) #p=0.0043. R2=59.93%


#Meta regression plot - linear
tiff(filename="SBP Meta Regression with Baseline SBP.tif", width=13.626, height=12, units='cm', res=300)
wi <- 1/sqrt(final$vi)
size <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))

plot(x=final$Baseline.BP, 
     y=final$yi, 
     xlab="Baseline SBP (mmHg)", 
     ylab="Mean Difference (mmHg)", 
     pch=19, 
     cex=size,
     ylim=c(-12, 2),
     xlim=c(100,150),
     bty='L') #datapoints sized by inverse variance

xi.new <- seq(100,160,length=100) #needed for regression line
pred <- predict(regBp, newmods = xi.new)
pred
lines(xi.new, pred$pred, lwd=2)
lines(xi.new, pred$ci.lb, lty="dashed", col="gray", lwd=2)
lines(xi.new, pred$ci.ub, lty="dashed", col="gray", lwd=2)
dev.off()
