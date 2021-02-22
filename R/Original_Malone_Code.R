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


library(raster);library(rgdal);library(tripack);library(SDMTools); library(manipulate);library(clhs);library(entropy)

setwd("Z:/Dropbox/2018/cLHC_samplingPAPER/gitHub/muddles/sampNumber") # set your working directory as appropriate


#data set (spatial data covariates)
df<- read.table("T1_covs.txt", header = T, sep = ",")
str(df)



####################################################################
# Data analysis for the population data

#Principal components of the population (the is for tests 1 and 2)
pca1 = prcomp(df[,3:6],scale=TRUE, center=TRUE)
scores_pca1 = as.data.frame(pca1$x)
screeplot(pca1) ## plot the variances explained by each component
biplot(pca1)
summary(pca1)

# retreive the loadings
pca1.load<- matrix(NA,ncol=4,nrow=4 )
for (i in 1:4){
  pca1.load[i,]<- as.matrix(t(pca1$rotation[i,]))}


#Quantiles of the population (this is for test 3)
# Number of bins
nb<- 25

#quantile matrix (of the covariate data)
q.mat<- matrix(NA, nrow=(nb+1), ncol= 4)
j=1
for (i in 3:ncol(df)){ #note the index start here
  #get a quantile matrix together of the covariates
  ran1<- max(df[,i]) - min(df[,i])
  step1<- ran1/nb 
  q.mat[,j]<- seq(min(df[,i]), to = max(df[,i]), by =step1)
  j<- j+1}
q.mat



#covariate data hypercube (this is for test 4)
## This takes a while to do so only do it once if you can 

cov.mat<- matrix(1, nrow=nb, ncol=4)
for (i in 1:nrow(df)){ # the number of pixels
  cntj<- 1 
  for (j in 3:ncol(df)){ #for each column
    dd<- df[i,j]  
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
cseq<- seq(10,500,10) # cLHC sample size
its<-10  # number internal iterations with each sample size number
mat.seq<- matrix(NA,ncol=8,nrow=length(cseq)) #empty matix for outputs


for (w in 1:length(cseq)){ # for every sample number configuration....
  s.size=cseq[w]  # sample size
  mat.f<- matrix(NA,ncol=8,nrow=its ) # placement for iteration outputs
  
  #internal loop
  for (j in 1:its){ #Note that this takes quite a while to run to completion
    repeat{
      ss <- clhs(df[,3:6], size = s.size, progress = T, iter = 100000) # Do a conditioned latin hypercube sample
      s.df<- df[ss,]
      if (sum(duplicated(s.df) | duplicated(s.df[nrow(s.df):1, ])[nrow(s.df):1]) < 2)
      {break}}
    
    #principal component of sample
    pca.s = prcomp(s.df[,3:6],scale=TRUE, center=TRUE)
    scores_pca1 = as.data.frame(pca.s$x)
    
    # plot the first 2 principal components and convex hull 
    rand.tr<-tri.mesh(scores_pca1[,1],scores_pca1[,2])
    rand.ch<-convex.hull(rand.tr, plot.it=F) #convex hull
    pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y)) # save the convext hull vertices
    #plot(scores_pca1[,1], scores_pca1[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])),ylim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])))
    #lines(c(rand.ch$x,rand.ch$x[1]), c(rand.ch$y,rand.ch$y[1]),col="red",lwd=1) # draw the convex hull(domain of prediction)
    
    #### First Test:
    #points in polygons routine
    # PCA prjection
    PCA_projection<- predict(pca.s, df[,3:6]) # Project population onto sample PC
    newScores = cbind(x=PCA_projection[,1],y=PCA_projection[,2]) # PC scores of projected population
    
    #plot the polygon and all points to be checked
    #plot(newScores,xlab="PCA 1", ylab="PCA 2", xlim=c(min(newScores[,1:2]), max(newScores[,1:2])),ylim=c(min(newScores[,1:2]), max(newScores[,1:2])),col='black', main='convex hull of ss')
    #polygon(pr_poly,col='#99999990')
    
    #create check which points fall within the polygon
    specMatch = pnt.in.poly(newScores,pr_poly)
    mat.f[j,7]<- sum(specMatch$pip)/nrow(specMatch)*100 # propertion of new spectra the fall within the convex hull
    #points(specMatch[which(specMatch$pip==0),1:2],pch='X', col='red')
    ##END points in polygons##
    
    #### Second Test:
    #similarity of the 2 matrices (PCA Similarity Factor; Krzanowski (1979))
    # retreive the loadings for the samples
    pca.s.load<- matrix(NA,ncol=4,nrow=4 )
    for (i in 1:4){
      pca.s.load[i,]<- as.matrix(t(pca.s$rotation[i,]))}
    
    # Perfrom the Krznowski 1979 calculation
    ps1<- pca1.load[,1:2]  
    ps2<- pca.s.load[,1:2]
    
    ps1.t<- t(ps1) #transpose
    ps2.t<- t(ps2) #transpose
    
    S<- ps1.t %*% ps2 %*% ps2.t %*% ps1
    mat.f[j,1]<-sum(diag(S))/2
    
    
    ## Third Test:
    #comparison of quantiles
    df.q2.1<-quantile(s.df$elevation, probs = seq(0, 1, 0.25),names = F, type = 7)
    df.q1.1<-quantile(df$elevation, probs = seq(0, 1, 0.25),names = F, type = 7)
    mat.f[j,2]<-sqrt((df.q1.1[1]-df.q2.1[1])^2 + (df.q1.1[2]-df.q2.1[2])^2 + (df.q1.1[3]-df.q2.1[3])^2 + (df.q1.1[4]-df.q2.1[4])^2 )
    
    df.q2.2<-quantile(s.df$yield, probs = seq(0, 1, 0.25),names = F, type = 7)
    df.q1.2<-quantile(df$yield, probs = seq(0, 1, 0.25),names = F, type = 7)
    mat.f[j,3]<-sqrt((df.q1.2[1]-df.q2.2[1])^2 + (df.q1.2[2]-df.q2.2[2])^2 + (df.q1.2[3]-df.q2.2[3])^2 + (df.q1.2[4]-df.q2.2[4])^2 )
    
    df.q2.3<-quantile(s.df$ECs, probs = seq(0, 1, 0.25),names = F, type = 7)
    df.q1.3<-quantile(df$ECs, probs = seq(0, 1, 0.25),names = F, type = 7)
    mat.f[j,4]<-sqrt((df.q1.3[1]-df.q2.3[1])^2 + (df.q1.3[2]-df.q2.3[2])^2 + (df.q1.3[3]-df.q2.3[3])^2 + (df.q1.3[4]-df.q2.3[4])^2 )
    
    df.q2.4<-quantile(s.df$ECd, probs = seq(0, 1, 0.25),names = F, type = 7)
    df.q1.4<-quantile(df$ECd, probs = seq(0, 1, 0.25),names = F, type = 7)
    mat.f[j,5]<-sqrt((df.q1.4[1]-df.q2.4[1])^2 + (df.q1.4[2]-df.q2.4[2])^2 + (df.q1.4[3]-df.q2.4[3])^2 + (df.q1.4[4]-df.q2.4[4])^2 )
    mat.f[j,6]<-mean(mat.f[j,2:5]) # take the mean distance
    
    
    ## Fourth test: Kullback-Leibler (KL) divergence
    ####Compare whole study area covariate space with the slected sample
    #sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
    h.mat<- matrix(1, nrow=nb, ncol=4)
    
    for (ii in 1:nrow(s.df)){ # the number of observations
      cntj<- 1 
      for (jj in 3:ncol(s.df)){ #for each column
        dd<- s.df[ii,jj]  
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
    klo<- mean(c(klo.1, klo.2,klo.3,klo.4))
    mat.f[j,8]<- klo  # value of 0 means no divergence
  } 
  
  
  #arrange outputs
  mat.seq[w,1]<-mean(mat.f[,6])
  mat.seq[w,2]<-sd(mat.f[,6])
  mat.seq[w,3]<-min(mat.f[,1])
  mat.seq[w,4]<-max(mat.f[,1])
  mat.seq[w,5]<-mean(mat.f[,7])
  mat.seq[w,6]<-sd(mat.f[,7])
  mat.seq[w,7]<-mean(mat.f[,8])
  mat.seq[w,8]<-sd(mat.f[,8])} ## END of LOOP

dat.seq<- as.data.frame(cbind(cseq,mat.seq))
names(dat.seq)<- c("samp_nos", "mean_dist","sd_dist", "min_S", "max_S", "mean_PIP","sd_PIP", "mean_KL","sd_KL")
##########################################################



#######################################################  
#plot some outputs
plot(cseq,mat.seq[,1], xlab="number of samples", ylab= "similarity between covariates (entire field) with covariates (sample)",main="Population and sample similarity")
plot(cseq,mat.seq[,2],xlab="number of samples", ylab= "standard deviation similarity between covariates (entire field) with covariates (sample)",main="Population and sample similarity (sd)")
plot(cseq,mat.seq[,3])
plot(cseq,mat.seq[,4])
plot(cseq,mat.seq[,5],xlab="number of samples", ylab= "percentage of total covariate variance of population account for in sample",main="Population and sample similarity")
plot(cseq,mat.seq[,6],xlab="number of samples", ylab= "standard deviation of percentage of total covariate variance of population account for in sample",main="Population and sample similarity")
plot(cseq,mat.seq[,7],xlab="number of samples", ylab= "KL divergence")
plot(cseq,mat.seq[,8],xlab="number of samples", ylab= "standard deviation of percentage of total covariate variance of population account for in sample",main="Population and sample similarity")
write.table(dat.seq, "Nav_datseq_clHC.txt", col.names=T, row.names=FALSE, sep=",")  # Save output to text file
##########################################################



##########################################################
# make an exponetial decay function (of the KL divergence)
x<- dat.seq$samp_nos
y = 1- (dat.seq$mean_PIP-min(dat.seq$mean_PIP))/(max(dat.seq$mean_PIP)-min(dat.seq$mean_PIP)) #PIP


#Parametise Exponential decay function
plot(x, y, xlab="sample number", ylab= "1 - PC similarity")          # Initial plot of the data
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
  b00=slider(0,1 , step=0.000001,initial= 0.01))

fit1 <- nls(y ~ k*exp(-b1*x) + b0, start = start)
summary(fit1) 
lines(x, fitted(fit1), col="red")


### Not used ###
#double exponential
#start <- list()
#manipulate(
#{
#  plot(x, y)
#  k <- kk; b0 <- b00; b1 <- b10; b2 <- b20
#  curve(k*(exp(-b1*x) + exp(-b2*x)) + b0, add=TRUE)
#  start <<- list(k=k, b0=b0, b1=b1, b2 = b2)
#},
#kk=slider(0, 5, step = 0.01,  initial = 1),
#b10=slider(0, 0.1, step = 0.000001, initial = 0.001),
#b20 = slider(0, 0.1, step = 0.00001, initial = 0.05),
#b00=slider(0,0.8 , step=0.00001,initial= 0.3))

#fit2 <- nls(y ~ k*(exp(-b1*x) + exp(-b2*x)) + b0, start = start, nls.control(maxiter = 1000))
#summary(fit2)
#lines(x, fitted(fit2), col="red")
#anova(fit1, fit2) 
### NOT USED ###
##############################################################################


#############################################################################
#Apply fit
xx<- seq(1, 500,1)
lines(xx, predict(fit1,list(x=xx)))

jj<- predict(fit1,list(x=xx))
normalized = 1- (jj-min(jj))/(max(jj)-min(jj))

x<- xx
y<- normalized

plot(x, y, xlab="sample number", ylab= "normalised PIP", type="l", lwd=2)          # Initial plot of the data

x1<- c(-1, 500); y1<- c(0.95, 0.95)
lines(x1,y1, lwd=2, col="red")

x2<- c(119, 119); y2<- c(0, 1)
lines(x2,y2, lwd=2, col="red")
#############################################################################

##END






