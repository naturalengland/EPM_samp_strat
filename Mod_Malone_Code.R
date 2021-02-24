# Script for determining the optimal number of samples to take using conditioned latin hypercube sampling.
# Given a suite of covariates this algorithm will assess an optimal number of samples based on a number of metrics.
# The metrics include:
# 1. Percentage of gird points within sample PC space
# 2. Principal component similarity
# 3. Quantile similarity
# 4. KL Divergence


######
# Note that in manuscript only the KL divergence is reported
# The follow script for each incremental sample size, the sample is repeat 10 times.
# We do this to look at the dispersion of the results resulting from different sample configurations of the same sample size
######


library(raster);library(rgdal);library(tripack);library(SDMTools); library(manipulate);library(clhs);library(entropy);library(ggplot2);library(ggpubr)

# setwd("F:/Projects/EPM/data/SWTest") # set your working directory as appropriate

setwd("../data/SWTest")

#generate Dataframe of covariates
####create stacks of Covariates & groups of Covariates####
CoV_SAR <- stack(lapply(list.files(path = "./AlignedRasters",
                                   pattern = paste0("*SAR_SW.tif$"),
                                   full.names = T), brick))

CoV_DTM <- stack(lapply(list.files(path = "./AlignedRasters", 
                                   pattern = paste0("*DTM_SW.tif$"), 
                                   full.names = T), raster))

CoV_Slope <- stack(lapply(list.files(path = "./AlignedRasters", 
                                   pattern = paste0("*Slope_SW.tif$"), 
                                   full.names = T), raster))

CoV_PSL <- stack(lapply(list.files(path = "./AlignedRasters", 
                                     pattern = paste0("*PeatySoils_SW.tif$"), 
                                     full.names = T), brick))

CoV_SupGeo <- stack(lapply(list.files(path = "./AlignedRasters", 
                                   pattern = paste0("*SupGeo_SW.tif$"), 
                                   full.names = T), raster))

CoV_All <- stack(CoV_SAR, CoV_DTM, CoV_Slope, CoV_PSL, CoV_SupGeo) #RasterStack of All of the CoVariates

CoV_sans_SAR <- stack(CoV_DTM, CoV_Slope, CoV_PSL, CoV_SupGeo) #RasterStack No SAR

CoV_sans_PSL <- stack(CoV_SAR, CoV_DTM, CoV_Slope, CoV_SupGeo) #RasterStack No PSL

CoV_Slope_Elev <- stack(CoV_DTM, CoV_Slope) #RasterStack Only Slope & elevation

#Check all covariates are loaded
names(CoV_All)
#creating a DataFrame of selected Covariates####
CoV_DF <- as.data.frame(rasterToPoints(CoV_All)) #change the Cov_x layer as required
str(CoV_DF)

####################################################################
# Data analysis for the population data

####USED for KL ####
#Quantiles of the population (this is for test 3)

# Number of bins
nb<- 25

#quantile matrix (of the covariate data)
q.mat<- matrix(NA, nrow=(nb+1), ncol= 9) #amended col to account for no. CoV.
j=1
for (i in 3:ncol(CoV_DF)){ #note the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(CoV_DF[,i]) - min(CoV_DF[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(CoV_DF[,i]), to = max(CoV_DF[,i]), by =step1)
  j<- j+1}
q.mat

#covariate data hypercube (this is for test 4)
## This takes a while to do so only do it once if you can < This really does take a long time.
cov.mat<- matrix(1, nrow=nb, ncol=9)
for (i in 1:nrow(CoV_DF)){ # the number of pixels
  cntj<- 1 
  for (j in 3:ncol(CoV_DF)){ #for each column #SHOULD WE ACTUALLY LEAVE IN THE X AND Y TO MAKE IT SPATIALLY EXPLICIT?
    dd<- CoV_DF[i,j]  
    for (k in 1:nb){  #for each quantile
      kl<- q.mat[k, cntj] 
      ku<- q.mat[k+1, cntj] 
      if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1} 
    }
    cntj<- cntj+1
  }
}

cov.mat
####################################################################





#######################################################################
#How many samples do we need?
#beginning of algorithm

#initial settings
cseq<- seq(50,500,50) # cLHC sample size, (beginning, end, step)
its.clhs <- 10 #number of iterations within clhs
its<-5  # number internal iterations with each sample size number
mat.seq<- matrix(NA,ncol=2,nrow=length(cseq)) #empty matrix for outputs
mat.f<- matrix(NA,ncol=2,nrow=its ) # placement for iteration outputs


for (w in 1:length(cseq)){ # for every sample number configuration....
  s.size=cseq[w]  # sample size
  print("Calculating KL divergence")
  
  #internal loop
  for (j in 1:its){ #Note that this takes quite a while to run to completion
    print(paste("evaluating sample of size ", s.size))
    print(paste("iteration", j, "of", its))
    repeat{
      start.rpt <- Sys.time()
      ss <- clhs(CoV_DF, size = s.size, progress = T, iter = its.clhs) # Do a conditioned latin hypercube sample
      s.CoV_DF<- CoV_DF[ss,]
      print(paste("time to calculate hypercube = ", lubridate::as.duration(Sys.time() - start.rpt)))
      if (sum(duplicated(s.CoV_DF) | duplicated(s.CoV_DF[nrow(s.CoV_DF):1, ])[nrow(s.CoV_DF):1]) < 2)
      {break}}
    

    ## Fourth test: Kullback-Leibler (KL) divergence####
    ####Compare whole study area covariate space with the slected sample
    #sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
    
    h.mat<- matrix(1, nrow=nb, ncol=9)
    
    for (ii in 1:nrow(s.CoV_DF)){ # the number of observations
      cntj<- 1 
      for (jj in 3:ncol(s.CoV_DF)){ #for each column
        dd<- s.CoV_DF[ii,jj]  
        for (kk in 1:nb){  #for each quantile
          kl<- q.mat[kk, cntj] 
          ku<- q.mat[kk+1, cntj] 
          if (dd >= kl & dd <= ku){h.mat[kk, cntj]<- h.mat[kk, cntj] + 1}
        }
        cntj<- cntj+1}}
    
    #h.mat 
    #Kullback-Leibler (KL) divergence
    klo.1<- KL.empirical(c(cov.mat[,1]), c(h.mat[,1])) #1
    klo.2<- KL.empirical(c(cov.mat[,2]), c(h.mat[,2])) #2
    klo.3<- KL.empirical(c(cov.mat[,3]), c(h.mat[,3])) #3
    klo.4<- KL.empirical(c(cov.mat[,4]), c(h.mat[,4])) #4
    klo.5<- KL.empirical(c(cov.mat[,5]), c(h.mat[,5])) #5
    klo.6<- KL.empirical(c(cov.mat[,6]), c(h.mat[,6])) #6
    klo.7<- KL.empirical(c(cov.mat[,7]), c(h.mat[,7])) #7
    klo.8<- KL.empirical(c(cov.mat[,8]), c(h.mat[,8])) #8
    klo.9<- KL.empirical(c(cov.mat[,9]), c(h.mat[,9])) #9

    klo<- mean(c(klo.1, klo.2,klo.3,klo.4,klo.5,klo.6,klo.7,klo.8,klo.9))
    mat.f[j,2]<- klo  # value of 0 means no divergence
    print(paste("KL divergence =", klo))
  } 
  
  
  #arrange outputs
  mat.seq[w,1]<-mean(mat.f[,2])
  mat.seq[w,2]<-sd(mat.f[,2])} ## END of LOOP

dat.seq<- as.data.frame(cbind(cseq,mat.seq))
names(dat.seq)<- c("samp_nos","mean_KL","sd_KL")
##########################################################



#######################################################  
#plot some outputs
plot(dat.seq[,1],dat.seq[,2],xlab="number of samples", ylab= "KL divergence")
plot(dat.seq[,1],dat.seq[,3],xlab="number of samples", ylab= "standard deviation of percentage of total covariate variance of population account for in sample",main="Population and sample similarity")
write.table(dat.seq, "Nav_datseq_clHC.txt", col.names=T, row.names=FALSE, sep=",")  # Save output to text file
##########################################################



##########################################################
#Individual Plots
#Plot 1 = Kl divergence by number of samples with line
x<- dat.seq$samp_nos
y <- dat.seq$mean_KL
start <- list()     # Initialize an empty list for the starting values


plot(x, y, xlab="sample number", ylab= "KL divergence")          # Initial plot of the data

#fit 1
manipulate(
  {
    plot(x, y)
    k <- kk; b0 <- b00; b1 <- b10
    curve(k*exp(-b1*x) + b0, add=TRUE)
    start <<- list(k=k, b0=b0, b1=b1)
  },
  kk=slider(0, 5, step = 0.01,  initial = 2),
  b10=slider(0, 1, step = 0.000001, initial = 0.01),
  b00=slider(0,1 , step=0.000001,initial= 0.01))

fit1 <- nls(y ~ k*exp(-b1*x) + b0, start = start)
summary(fit1) 
lines(x, fitted(fit1), col="red")

#DF OF PREDICTED VALUES FOR SAMPLE SIZE
Pred <- data.frame(2:(max(cseq)+1))
Pred$samp_nos <- seq(1, max(cseq), 1)
Pred$Pred_KL <- predict(fit1,list(x=Pred$samp_nos))
Pred [1] <- NULL

#Plot
KL_Plot <- ggplot() #create plot frame

KL_Plot <- KL_Plot + #add the data points
  geom_point(aes(x = as.numeric(dat.seq$samp_nos), y = as.numeric(dat.seq$mean_KL))) +
  xlab("Sample Number") +
  ylab("KL divergence")

KL_Plot <- KL_Plot + geom_line(aes(x = Pred$samp_nos, y = Pred$Pred_KL), colour = "blue") +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.05)) +
  scale_x_continuous(limits = c(0, max(cseq)), breaks = seq(0, max(cseq), (max(cseq)/5))) +
  theme(legend.position = "none") #add the fitted line

KL_Plot

#########################################

#Plot 2 cdf of 1-KL (not normalised)
x <- dat.seq$samp_nos
y2 <- 1 - dat.seq$mean_KL

plot(x, y2, xlab="sample number", ylab= "CDF 1-KL")         # Initial plot of the data

#fit 2
manipulate(
  {
    plot(x, y2)
    k2 <- kk2; b02 <- b002; b12 <- b102
    curve(k2*exp(-b12*x) + b02, add=TRUE)
    start <<- list(k2=k2, b02=b02, b12=b12)
  },
  kk2=slider(0, 5, step = 0.01,  initial = 2),
  b102=slider(0, 1, step = 0.000001, initial = 0.01),
  b002=slider(0,1 , step=0.000001,initial= 0.01))

fit2 <- nls(y2 ~ k2*exp(-b12*x) + b02, start = start)
summary(fit2) 
lines(x, fitted(fit2), col="red")

#DF OF PREDICTED CDF VALUES FOR SAMPLE SIZE
PredCDF <- data.frame(2:(max(cseq)+1))
PredCDF$samp_nos <- seq(1, max(cseq), 1)
PredCDF$Pred_CDF <- predict(fit2,list(x=PredCDF$samp_nos))
PredCDF [1] <- NULL

ConfAchieved <- as.numeric(which(PredCDF$Pred_CDF >= 0.95) [1])

#Plot
CDF_Plot <- ggplot() #create plot frame

CDF_Plot <- CDF_Plot + #add the data points
  geom_point(aes(x = as.numeric(dat.seq$samp_nos), y = 1- as.numeric(dat.seq$mean_KL))) +
  xlab("Sample Number") +
  ylab("CDF of 1- KL divergence")

CDF_Plot <- CDF_Plot + geom_line(aes(x = PredCDF$samp_nos, y = PredCDF$Pred_CDF), colour = "blue") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, max(cseq)), breaks = seq(0, max(cseq), (max(cseq)/5))) +
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0.95, linetype = "dashed", alpha = 0.4, colour = "red") +
  geom_vline(xintercept = ConfAchieved, linetype = "dashed", alpha = 0.4, colour = "red") +
  geom_text(aes(x=ConfAchieved, label= paste0(ConfAchieved,"\n"), y=0.2), colour="red", angle=90)


CDF_Plot

##Combine plots to one figure####
All_Cov_Fig <- ggarrange(KL_Plot, CDF_Plot,
                         labels = c("A", "B"),
                         ncol = 1, nrow = 2)

All_Cov_Fig
####################################
##END