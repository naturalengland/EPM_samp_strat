######
# Based on Malone et al. code.
# expanded upon to accept both numerical and factor covariates
# always keep numerical covariates at the start of the dataframe followed by factor covariates
# changes to code necessary to allow CoV classes to be in any order.
######


library(raster);library(rgdal);library(tripack);#library(SDMTools); 
library(manipulate);library(clhs);library(entropy);library(ggplot2);
library(ggpubr);library(openxlsx); library(plyr)


#setwd("F:/Projects/EPM/data/AlignedNationWide") # set your working directory as appropriate

setwd("../data")

#Automation, for a more powerful computer/ cloud?
for (BGZ in c(1:14)) { #this will be 14 when we run this properly
  BGZ.St <- Sys.time()
  print(paste0("Starting Analysis of BGZ ",BGZ))
  #BGZ <- 2 #usually comment this out!

  #generate Dataframe of covariates
  ####create stacks of Covariates & groups of Covariates####
  CoV_SAR <- stack(lapply(list.files(path = "./AlignedNationWide",
                                     pattern = paste0("SAR_", BGZ,".tif$"),
                                     full.names = T, recursive = T), brick))
  
  CoV_DTM <- stack(lapply(list.files(path = "./AlignedNationWide",
                                     pattern = paste0("DTM_", BGZ,".tif$"), 
                                     full.names = T, recursive = T), raster))
  
  CoV_Slope <- stack(lapply(list.files(path = "./AlignedNationWide",
                                       pattern = paste0("Slope_", BGZ,".tif$"), 
                                       full.names = T, recursive = T), raster))
  
  CoV_PSL <- stack(lapply(list.files(path = "./AlignedNationWide",
                                     pattern = paste0("PSL_", BGZ,".tif$"),
                                     full.names = T, recursive = T), brick))
  
  CoV_SupGeo <- stack(lapply(list.files(path = "./AlignedNationWide",
                                        pattern = paste0("SupGeo_", BGZ,".tif$"),
                                        full.names = T, recursive = T), raster))
  
  CoV_BedGeo <- stack(lapply(list.files(path = "./AlignedNationWide",
                                        pattern = paste0("BedGeo_", BGZ,".tif$"), 
                                        full.names = T, recursive = T), raster))
  
  CoV_All <- stack(CoV_SAR, CoV_DTM, CoV_Slope, CoV_PSL, CoV_SupGeo, CoV_BedGeo) #RasterStack of All of the CoVariates
  
#########################################
#Check all covariates are loaded
names(CoV_All)

#creating a DataFrame of Covariates####
CoV_DF <- as.data.frame(rasterToPoints(CoV_All)) #change the Cov_x layer as required
CoV_DF[, grep("SupGeo", colnames(CoV_DF))] <- as.factor(CoV_DF[, grep("SupGeo", colnames(CoV_DF))])
CoV_DF[, grep("Bed", colnames(CoV_DF))] <- as.factor(CoV_DF[, grep("Bed", colnames(CoV_DF))])

#check for and remove columns that have only 1 level e.g all 0 (no coverage) in one of the PSLs
no_uniques <- function(vec) length(unique(na.omit(vec))) > 1

CoV_DF <- CoV_DF[,sapply(CoV_DF,no_uniques)]
str(CoV_DF)

# #select everything but SAR
# CoV_sans_SAR <- CoV_DF[, -grep("SAR_...$", names(CoV_DF))]
# 
# #select everything but PSL
# CoV_sans_PSL <- CoV_DF[, -grep("PSL_...$", names(CoV_DF))]
# 
# #select only slope & DTM
# CoV_Slope_Elev <- CoV_DF[, grep("^x$|^y$|Slope|DTM", names(CoV_DF))]

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
  ran1<- max(na.omit(CoV_DF[,i])) - min(na.omit(CoV_DF[,i]))
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(na.omit(CoV_DF[,i])), to = max(na.omit(CoV_DF[,i])), by =step1)
  j<- j+1}

colnames(q.mat)<- c(num_CoV_list)

############################################
#creating a matrix for the factor CoVs, no longer used in calcs but produced for output check
fac.mat<- matrix(NA, nrow= max_levels, ncol= no_factors)
j=1
for (i in 1:length(fac_levels)){ #return the col. number of each factor variable
  fac.mat[1:fac_levels[[j]],j] <- as.numeric(levels
                                             (subset(CoV_DF, select = fac_CoV_cols)[,i])) # print the levels of each factor into the matrix
  j<- j+1}

colnames(fac.mat)<- c(fac_CoV_list)

############################################

CoVntiles <- CoV_DF[num_CoV_cols]

CoV_counts <- data.frame(matrix(NA, nrow = nb, ncol = 0))

CoV_Col <- 1
Num_CoV <- 1

for (Num_CoV in 1:no_cont){
#compute a histogram for the breaks as specified in q.mat
   Cov_hist <- hist(as.matrix(CoVntiles[,Num_CoV]), 
                 breaks= seq(min(q.mat[,CoV_Col]),
                             max(q.mat[,CoV_Col]),
                 by=(as.numeric(q.mat[2,Num_CoV])
                     - as.numeric(q.mat[1,Num_CoV]))), plot = FALSE)

   #add this to the dataframe of counts.
   CoV_counts[ , ncol(CoV_counts) + 1] <- Cov_hist$counts                              # Append new column
   colnames(CoV_counts)[ncol(CoV_counts)] <- paste0(names(CoVntiles[CoV_Col]))         # Rename column name

CoV_Col <- CoV_Col + 1
}

cov.mat <- as.matrix(CoV_counts+1)
##########################################

Facntiles <- CoV_DF[fac_CoV_cols]

#get the levels of each factor and create a matrix list from them
Get_levels.m <- matrix(sapply(CoV_DF[fac_CoV_cols], levels), ncol = no_factors, byrow = TRUE)

# #list of levels in the factor
# Get_levels.m[[1]]

#empty vector for counts
Fac_counts.v <- c()

for (f.i in 1:no_factors){
  #make a vector to add the count each time to
  Fac_counts.v <- c(Fac_counts.v, 
                      list(as.data.frame(table(
                      factor(Facntiles[,f.i],
                      levels = Get_levels.m[[f.i]])))[,2]))
}


# #For Output
# fac.cov.mat<- matrix(NA, nrow= max_levels, ncol= no_factors)
# fac.cov.mat.1 <- cbind(unlist(Fac_counts.v[1]))
# fac.cov.mat.2 <- cbind(unlist(Fac.Counts[2]))

#HACKY FOR TEST
CoVGroups <- list(CoV_DF)
Group <- 1

#made to work for code review this is just maintaining the matrices created above.
q.mat.gr <- q.mat
fac.mat.gr <- fac.mat
cov.mat.gr <- cov.mat
#######################################################################
#How many samples do we need?
#beginning of algorithm

#initial settings
cseq<- seq(100, 350, 50) # cLHC sample size, (beginning, end, step)
its.clhs <- 100 #number of iterations within clhs
its<-5  # number internal iterations with each sample size number
mat.seq<- matrix(NA,ncol=2,nrow=length(cseq)) #empty matrix for outputs
mat.f<- matrix(NA,ncol=2,nrow=its ) # placement for iteration outputs


for (w in 1:length(cseq)){ # for every sample number configuration....
  s.size=cseq[w]  # sample size
  
  #internal loop
  for (j in 1:its){ #Note that this takes quite a while to run to completion
    print(paste("evaluating sample of size ", s.size))
    print(paste("iteration", j, "of", its))
    repeat{
      start.rpt <- Sys.time()
      ss <- clhs(CoVGroups[[Group]][,3:no_CoV], size = s.size, progress = T, iter = its.clhs) # Do a conditioned latin hypercube sample
      s.CoV_DF<- CoVGroups[[Group]][ss,] #select the row numbers output by clhs from the cov.
      print(paste("time to calculate hypercube = ", lubridate::as.duration(Sys.time() - start.rpt)))
      if (sum(duplicated(s.CoV_DF) | duplicated(s.CoV_DF[nrow(s.CoV_DF):1, ])[nrow(s.CoV_DF):1]) < 2)
      {break}}
    
    ############################################
    ## Fourth test: Kullback-Leibler (KL) divergence####
    ####Compare whole study area covariate space with the slected sample
    #sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)

    ############################################
    s.CoVntiles <- s.CoV_DF[num_CoV_cols]
    
    #empty dataframe for counts
    h.mat <- data.frame(matrix(NA, nrow = nb, ncol = 0))
    
    #resetting these vars just in case/ for testing
    CoV_Col <- 1
    Num_CoV <- 1
    
    #for each continuous covariate make a histogram with the bounds as set by q.mat
    for (Num_CoV in 1:no_cont){
      #group into quantiles
      Cov_hist <- hist(as.matrix(s.CoVntiles[,Num_CoV]), 
                       breaks= seq(min(q.mat[,CoV_Col]),
                                   max(q.mat[,CoV_Col]),
                                   by=(as.numeric(q.mat[2,Num_CoV])
                                       - as.numeric(q.mat[1,Num_CoV]))), plot = FALSE)
      
      
      #take the counts from the histogram and add them to the dataframe created earlier
      h.mat[ , ncol(h.mat) + 1] <- Cov_hist$counts                              # Append new column
      colnames(h.mat)[ncol(h.mat)] <- paste0(names(s.CoVntiles[CoV_Col]))         # Rename column name
      
      CoV_Col <- CoV_Col + 1
    }
    
    #transform to matrix
    h.mat <- as.matrix(h.mat+1)
    
    
    ############################################
    #create a dataframe of all of the factor covariate values in the sample
    h.Facntiles <- s.CoV_DF[fac_CoV_cols]
    
    #empty vector for counts
    h_Fac_counts.v <- c()
    
    for (h.f.x in 1:no_factors){
      #count the occurance of each factor level and add as list to vector
      h_Fac_counts.v <- c(h_Fac_counts.v, 
                          list(as.data.frame(table(
                          factor(h.Facntiles[,h.f.x],
                          levels = Get_levels.m[[h.f.x]])))[,2]))
    } 
    
    
    ############################################
    
    #Kullback-Leibler (KL) divergence
    #
    klo.v <- c() #empty vector for covar results
    for (iiii in num_CoV_cols){
      n.x <- iiii -2
      klo.v <- c(klo.v, #bind following result to vector
                 KL.empirical(cov.mat.gr[,n.x], h.mat[,n.x])) #loop for numerical  
    }
    
    for (jjjj in fac_CoV_cols) {
      f.x <- (jjjj-no_cont)-2
      klo.v <- c(klo.v, #bind following result to vector
                 KL.empirical(c(unlist(Fac_counts.v[f.x])+1),
                                c(unlist(h_Fac_counts.v[f.x])+1)))
                   
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
x <- dat.seq$samp_nos
y <- dat.seq$mean_KL

x.samp.rng <- seq(min(cseq),max(cseq),1)

plot(x, y, xlab="sample number", ylab= "KL divergence")          # Initial plot of the data
#########################################
# fitlog1 <- lm(y~log(x))
# summary(fitlog1)
# 
# Pred1 <- predict(fitlog1, newdata = list(x = x.samp.rng) ,interval ="confidence")
# 
# #DF OF PREDICTED VALUES FOR SAMPLE SIZE
# Pred1.df <- data.frame(Pred1)
# Pred1.df$samp_nos <- seq(min(cseq),max(cseq),1)
# 
# #Plot
# KL_Plot <- ggplot() #create plot frame
# 
# KL_Plot <- KL_Plot + #add the data points
#   geom_point(aes(x = as.numeric(dat.seq$samp_nos), y = as.numeric(dat.seq$mean_KL))) +
#   xlab("Sample Number") +
#   ylab("KL divergence")
# 
# KL_Plot <- KL_Plot + geom_line(aes(x = x.samp.rng, y = Pred1.df$fit), 
#                                colour = "black") +
#                      geom_line(aes(x = Pred1.df$samp_nos, y = Pred1.df$upr), 
#                                colour = "red", 
#                                linetype = "longdash") +
#                      geom_line(aes(x = Pred1.df$samp_nos, y = Pred1.df$lwr), 
#                                colour = "green", 
#                                linetype = "longdash") +
#                      scale_y_continuous(limits = c(0, 0.4), 
#                                          breaks = seq(0, 0.4, 0.05))
#   
#   
# KL_Plot

#########################################
start <- list()     # Initialize an empty list for the starting values
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
#########################################
# fitlog2 <- lm(y2~log(x))
# summary(fitlog2)
# 
# Pred2 <- predict(fitlog2, newdata = list(x = x.samp.rng) ,interval ="confidence")
# 
# #DF OF PREDICTED VALUES FOR SAMPLE SIZE
# Pred2.df <- data.frame(Pred2)
# Pred2.df$samp_nos <- seq(min(cseq),max(cseq),1)
# 
# ConfAchieved <- Pred2.df[as.numeric(min(which(Pred2.df$fit >= .95))),4]
# 
# 
# #Plot
# CDF_Plot <- ggplot() #create plot frame
# 
# CDF_Plot <- CDF_Plot + #add the data points
#   geom_point(aes(x = as.numeric(dat.seq$samp_nos), y = 1- as.numeric(dat.seq$mean_KL))) +
#   xlab("Sample Number") +
#   ylab("CDF of 1- KL divergence")
# 
# CDF_Plot <- CDF_Plot + geom_line(aes(x = x.samp.rng, y = Pred2.df$fit), 
#                                  colour = "black") +
#                        geom_line(aes(x = Pred2.df$samp_nos, y = Pred2.df$upr), 
#                                  colour = "red", 
#                                  linetype = "longdash") +
#                        geom_line(aes(x = Pred2.df$samp_nos, y = Pred2.df$lwr), 
#                                  colour = "green", 
#                                  linetype = "longdash") + 
#   geom_hline(yintercept = 0.95, linetype = "dashed", alpha = 0.4, colour = "blue") +
#   geom_vline(xintercept = ConfAchieved, linetype = "dashed", alpha = 0.4, colour = "blue") +
#   geom_text(aes(x=ConfAchieved, label= paste0(ConfAchieved,"\n"), y=0.2), colour="blue", angle=90)
# 
# 
# CDF_Plot
#########################################
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

wb <- createWorkbook()

addWorksheet(wb, "q.mat.gr")
writeData(wb, sheet = "q.mat.gr", rowNames = TRUE, colNames = TRUE, x = q.mat.gr)

addWorksheet(wb, "fac.mat.gr")
writeData(wb, sheet = "fac.mat.gr", rowNames = TRUE, colNames = TRUE, x = fac.mat.gr)

addWorksheet(wb, "cov.mat.gr")
writeData(wb, sheet = "cov.mat.gr", rowNames = TRUE, colNames = TRUE, x = cov.mat.gr)

# addWorksheet(wb, "fac.cov.mat.1")
# writeData(wb, sheet = "fac.cov.mat.1", rowNames = TRUE, colNames = TRUE, x = fac.cov.mat.1)
# 
# addWorksheet(wb, "fac.cov.mat.2")
# writeData(wb, sheet = "fac.cov.mat.2", rowNames = TRUE, colNames = TRUE, x = fac.cov.mat.2)

addWorksheet(wb, "h.mat")
writeData(wb, sheet = "h.mat", rowNames = TRUE, colNames = TRUE, x = h.mat)

# addWorksheet(wb, "h.fac.cov.mat.1")
# writeData(wb, sheet = "h.fac.cov.mat.1", rowNames = TRUE, colNames = TRUE, x = h.fac.cov.mat.1)
# 
# addWorksheet(wb, "h.fac.cov.mat.2")
# writeData(wb, sheet = "h.fac.cov.mat.2", rowNames = TRUE, colNames = TRUE, x = h.fac.cov.mat.2)

addWorksheet(wb, "dat.seq")
writeData(wb, sheet = "dat.seq", rowNames = TRUE, colNames = TRUE, x = dat.seq)

addWorksheet(wb, "ConfAchieved")
writeData(wb, sheet = "ConfAchieved", x = ConfAchieved)

saveWorkbook(wb, file=paste0("./ModelOutputs/",BGZ,"/BGZ_",BGZ,".xlsx"))
#./ModelOutputs/",BGZ,"/Group ",Group,"/BGZ_",BGZ,"_",Group,".xlsx"))

# 1. Openfile
png(file=paste0("./ModelOutputs/",BGZ,"/BGZ_",BGZ,".png")
)

# 2. Create the plot
plot(All_Cov_Fig)

# 3. Close the file
dev.off()


# 1. Open file
pdf(file=paste0("./ModelOutputs/",BGZ,"/BGZ_",BGZ,".pdf"), 
    width = 10, 
    height = 8
)
# 2. Create the plot
plot(All_Cov_Fig)

# 3. Close the file
dev.off()

#}

print(paste("Time to Complete BGZ",BGZ," = ", lubridate::as.duration(Sys.time() - BGZ.St)))

#clear the variables so that there is sufficient resources to remake for next run
rm(list=ls())

}
