
######
# Based on Malone et al. code.
# expanded upon to accept both numerical and factor covariates
# always keep numerical covariates at the start of the dataframe followed by factor covariates
# changes to code necessary to allow CoV classes to be in any order.
######


library(raster);library(rgdal);library(tripack);library(SDMTools); library(manipulate);library(clhs);library(entropy);library(ggplot2);library(ggpubr)

setwd("E:/Natural England/Peatland/EPM_samp_strat/") # set your working directory as appropriate

#setwd("")

#generate Dataframe of covariates
####create stacks of Covariates & groups of Covariates####
CoV_SAR <- stack(lapply(list.files(path = "./data",
                                   pattern = paste0("*SAR.tif$"),
                                   full.names = T), brick))

CoV_DTM <- stack(lapply(list.files(path = "./data",
                                   pattern = paste0("*DTM.tif$"), 
                                   full.names = T), raster))

CoV_Slope <- stack(lapply(list.files(path = "./data",
                                     pattern = paste0("*Slope.tif$"), 
                                     full.names = T), raster))

CoV_PSL <- stack(lapply(list.files(path = "./data",
                                   pattern = paste0("*PSL.tif$"), 
                                   full.names = T), brick))

CoV_SupGeo <- stack(lapply(list.files(path = "./data",
                                      pattern = paste0("*SupGeo.tif$"), 
                                      full.names = T), raster))

CoV_BedGeo <- stack(lapply(list.files(path = "./data",
                                      pattern = paste0("*Bed_Geo.tif$"), 
                                      full.names = T), raster))

CoV_All <- stack(CoV_SAR, CoV_DTM, CoV_Slope, CoV_PSL, CoV_SupGeo, CoV_BedGeo) #RasterStack of All of the CoVariates

#########################################
#Check all covariates are loaded
names(CoV_All)

#creating a DataFrame of selected Covariates####
CoV_DF <- as.data.frame(rasterToPoints(CoV_All)) #change the Cov_x layer as required
CoV_DF[, grep("SupGeo", colnames(CoV_DF))] <- as.factor(CoV_DF[, grep("SupGeo", colnames(CoV_DF))])
CoV_DF[, grep("Bed", colnames(CoV_DF))] <- as.factor(CoV_DF[, grep("Bed", colnames(CoV_DF))])
str(CoV_DF)

CoV_sans_SAR <- subset(CoV_DF, select =c(1,2,6:ncol(CoV_DF)))
str(CoV_sans_SAR)

CoV_sans_PSL <- subset(CoV_DF, select =c(1:7,11,12)) #RasterStack No PSL
str(CoV_sans_PSL)

CoV_Slope_Elev <- subset(CoV_DF, select =c(1,2,6,7))
str(CoV_Slope_Elev)

####################################################################
# Data analysis for the population data
#create a list of the factor covs
fac_CoV_list <- as.list(names(CoV_DF)[sapply(CoV_DF, is.factor)])

#return the position of those columns within the DF
fac_CoV_cols <- unlist(lapply(fac_CoV_list, function(x) which(colnames(CoV_DF) == x)))

#create a list of the numerical covs
num_CoV_list <- as.list(names(CoV_DF)[sapply(CoV_DF, is.numeric)])[-c(1,2)]

#return the position of those columns within the DF
num_CoV_cols <- unlist(lapply(num_CoV_list, function(x) which(colnames(CoV_DF) == x)))

#useful variables, each self explanatory.
no_CoV <- length(names(CoV_DF))-2 #account for removal of x y.
no_cont <- length(names(CoV_DF)[sapply(CoV_DF, is.numeric)])-2 #account for removal of x y.
no_factors <- length(names(CoV_DF)[sapply(CoV_DF, is.factor)])
max_levels <- max(sapply(seq(1, ncol(CoV_DF)), 
                         function(x) {length(levels(CoV_DF[,x]))}))
fac_levels <- lapply(fac_CoV_list, function(x) length(levels(CoV_DF[,x])))

# Number of bins
nb <- 25
############################################
##creating a quantile matrix for continuous vars
q.mat<- matrix(NA, nrow=(nb+1), ncol= no_cont) 
j=1
for (i in num_CoV_cols){ #note the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(CoV_DF[,i]) - min(CoV_DF[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(CoV_DF[,i]), to = max(CoV_DF[,i]), by =step1)
  j<- j+1}
q.mat

############################################
#creating a matrix for the factor CoVs
fac.mat<- matrix(NA, nrow= max_levels, ncol= no_factors)
j=1
for (i in 1:length(fac_levels)){ #return the col. number of each factor variable
  fac.mat[1:fac_levels[[j]],j] <- as.numeric(levels
                                             (subset(CoV_DF, select = fac_CoV_cols)[,i])) # print the levels of each factor into the matrix
  j<- j+1}

fac.mat

fac_levels[[1]]

############################################
#covariate data hypercube (this is for test 4)
## This takes a while to do so only do it once if you can < This really does take a long time.
cov.mat<- matrix(1, nrow=nb, ncol=no_cont)
for (i in 1:nrow(CoV_DF)){ # the number of pixels 
  cntj<- 1 
  for (j in num_CoV_cols){ #for each column
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

############################################
#this produces a new covariate matrix for the factor covariates.
fac.cov.mat<- matrix(1, 
                     nrow = max_levels,
                     ncol = no_factors)
for (i in 1:nrow(CoV_DF)){ # the number of pixels 
  cntj<- 1 
  for (j in fac_CoV_cols){ #for each column that contains a factor covariate
    dd<- CoV_DF[i,j]  
    for (k in 1:length(levels(CoV_DF[,j]))){  #for each variable, the number of buckets (levels)
      kl<- fac.mat[k, cntj] 
      if (dd == kl){fac.cov.mat[k, cntj]<- fac.cov.mat[k, cntj] + 1} #add 1 to the count if the value is equal to the bucket (level)
    }
    cntj<- cntj+1
  }
}
fac.cov.mat
############################################
####################################################################





#######################################################################
#How many samples do we need?
#beginning of algorithm

#initial settings
cseq<- seq(100,500,10) # cLHC sample size, (beginning, end, step)
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
      s.CoV_DF<- CoV_DF[ss,] #select the row numbers output by clhs from the cov.
      print(paste("time to calculate hypercube = ", lubridate::as.duration(Sys.time() - start.rpt)))
      if (sum(duplicated(s.CoV_DF) | duplicated(s.CoV_DF[nrow(s.CoV_DF):1, ])[nrow(s.CoV_DF):1]) < 2)
      {break}}
    
    ############################################
    ## Fourth test: Kullback-Leibler (KL) divergence####
    ####Compare whole study area covariate space with the slected sample
    #sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
    
    h.mat<- matrix(1, nrow=nb, ncol=no_cont)
    
    for (ii in 1:nrow(s.CoV_DF)){ # the number of pixels in sample
      cntj<- 1 
      for (jj in num_CoV_cols){ #for each column of cont (numeric) vars.
        dd<- s.CoV_DF[ii,jj]
        for (kk in 1:nb){  #for each quantile
          kl<- q.mat[kk, cntj] 
          ku<- q.mat[kk+1, cntj] 
          if (dd >= kl & dd <= ku){h.mat[kk, cntj]<- h.mat[kk, cntj] + 1} 
        }
        cntj<- cntj+1
      }
    }
    
    h.mat
    
    ############################################
    #this produces a new covariate matrix for the factor covariates.
    h.fac.mat<- matrix(1, 
                       nrow = max_levels,
                       ncol = no_factors)
    for (iii in 1:nrow(s.CoV_DF)){ # the number of pixels in sample
      cntj<- 1 
      for (jjj in fac_CoV_cols){ #for each column
        dd<- s.CoV_DF[iii,jjj]  #get the value of each in covariate for each sample
        for (kkk in 1:length(levels(s.CoV_DF[,jjj]))){  #for each variable the number of buckets (levels)
          kl<- fac.mat[kkk, cntj] 
          if (dd == kl){h.fac.mat[kkk, cntj]<- h.fac.mat[kkk, cntj] + 1} 
        }
        cntj<- cntj+1
      }
    }
    h.fac.mat
    
    ############################################ TO HERE>
    
    #Kullback-Leibler (KL) divergence
    #
    klo.v <- c() #empty vector for covar results
    for (iiii in num_CoV_cols){
      x <- iiii -2
      klo.v <- c(klo.v, #bind following result to vector
                 KL.empirical(c(cov.mat[,x]), c(h.mat[,x]))) #loop for numerical  
    }
    
    for (jjjj in fac_CoV_cols) {
      x <- (jjjj-no_cont)-2
      klo.v <- c(klo.v, #bind following result to vector
                 KL.empirical(c(fac.cov.mat[1:fac_levels[[x]],x]),
                              c(h.fac.mat[1:fac_levels[[x]],x]))) #loop for factors, comparing only the levels of data not additional rows in matrix.
    }
    
    
    klo<- mean(klo.v)
    mat.f[j,2]<- klo  # value of 0 means no divergence
    print(paste("KL divergence =", klo))
  }
  
  
  #arrange outputs
  mat.seq[w,1]<-mean(mat.f[,2])
  mat.seq[w,2]<-sd(mat.f[,2])} ## END of LOOP

dat.seq<- as.data.frame(cbind(cseq,mat.seq)) #writes the number of samples next to the outputs
names(dat.seq)<- c("samp_nos","mean_KL","sd_KL")
dat.seq
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
  b00=slider(0, 1, step=0.000001,initial= 0.01))

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