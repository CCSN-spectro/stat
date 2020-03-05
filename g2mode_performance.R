source("specPdgrm.R")
source("data_generator.R")
source("functions.R")
library(signal)

##################
### parameters ###
##################
fs=4096
duration=1
ampl=10
fcut=15

######################
### prepare signal ###
######################
s20 = read.table("../simulation_data/s20-gw_10kpc16384.dat"); # V1 time, V2 signal

# signal sampled at 16384 Hz. Resampling at fs
signaly=resample(s20$V2,1/16384,1/fs)
# decimate time vector
signalx=seq(1,fs*duration,by=1)
for (i in 1:fs*duration) {
  signalx[i]=mean(s20$V1[((i-1)*4+1):((i-1)*4+4)])
}

wvf.df = data.frame("V1"=signalx,"V2"=signaly)

# determine the time interval where signal is non null
s20_0 = s20[,2]!=0;  # identifying values == 0
s20_0 = s20[s20_0,]; # discarding values

timeGmode = s20_0[c(1,length(s20_0[,1])), 1];#time interval for signal different from 0

###################
### True ratios ###
###################

# For g2 modes, ratio (x variable in TF19) is x = Mpns / Rpns^2 
true_data = read.table("../simulation_data/s20_data_g2.dat", sep = ",", comment.char = "#",header = TRUE); # true values
colnames(true_data) = c ("time","x");

############################
### FREQUENTIST ANALYSIS ###
############################

#############
### Model ###
#############

#fits_data = read.table("inputs/fits_data.dat", sep = ",");# data to generate model
fits_data = read.table("../simulation_data/A-A_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("../simulation_data/CoCo_fits_data_g2.dat", sep = ",");# data to generate model
colnames(fits_data) = c("r", "f");


# Linear model
mod = lm(r~ 0 + f + I(f^3), data = fits_data);# final model
summary(mod);

# 
mu1 = "~ f - 1";
mu2 = "~ f + I(f^2) - 1";
mu3 = "~ f + I(f^2) + I(f^3) - 1";
mu4 = "~ f + I(f^3) - 1";
s1 = "~ f - 1"; 
s2 = "~ f + I(f^2) - 1";

Xm  = model.matrix(eval(parse(text=eval(parse(text=mu2)))), fits_data);
Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = FALSE);

# g-mode estimate: starting from the right, left or median
#gmode = c("left", "right", "median");
gmode = c("left");
#gmode = c("right"); 

###########
### Run ###
###########

# loop over N generation of noisy data and add signal

N=100
dist_nb=24
result<-matrix(nrow=N*dist_nb,ncol=12)

for (j in 1:dist_nb){
  dist = 1+(j-1)*.25
  print(c("distance:",dist))
  for (i in 1:N){
    d = data_generator(fs, duration, wvf.df, ampl=10/dist, filtering = "whitening_bis", fcut, actPlot=FALSE);
    noisydata = data.frame("V1"=d$t,"V2"=d$y);
    out = covpbb(noisydata, mod=mod, l=200, p=90, fs=fs,
              um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(-Inf, 500),
              um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(500, 1700),
              gmode = gmode,
              thruth_data=true_data, actPlot=FALSE,
              limFreq = c(1000));

    out1 = covpbb(noisydata, mod=fit, l=200, p=90, fs=fs,
               um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(-Inf, Inf),
               um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(500, 1700),
               gmode = gmode,
               thruth_data=true_data, actPlot=FALSE,
               limFreq = c(1000));
    
    result[i+(j-1)*N,1]=out$covpbb[1,1]
    result[i+(j-1)*N,2]=dist

    result[i+(j-1)*N,3]=out$covpbb[1,2]
    result[i+(j-1)*N,4]=out$covpbb[1,3]
    result[i+(j-1)*N,5]=out$residual[1,2]
    result[i+(j-1)*N,6]=out$residual[1,3]
    result[i+(j-1)*N,7]=out$chi2[1,2]
  
    result[i+(j-1)*N,8]=out1$covpbb[1,2]
    result[i+(j-1)*N,9]=out1$covpbb[1,3]
    result[i+(j-1)*N,10]=out1$residual[1,2]
    result[i+(j-1)*N,11]=out1$residual[1,3]
    result[i+(j-1)*N,12]=out1$chi2[1,2]
  
  }
}
write.table(result, file="results_fcut1000Hz_AA.txt", sep=" ", row.names=FALSE, col.names=FALSE)

#########################
### BAYESIAN ANALYSIS ###
#########################
  
# postSamples = read.table("postsamples.txt", header = T); # importing posterios samples
# 
# covpbbBayes(noisydata, postSamples, l=200, p=90, fs=fs, movGmode = 11, 
#             um = 10, dm = 10, movBand = 5, timeGmode = NULL, 
#             Obergaunlinguer_data, actPlot =TRUE, fmean = mean(fits_data$f))
