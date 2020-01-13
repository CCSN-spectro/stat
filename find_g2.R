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
#s20 = read.table("../simulation_data/s20-gw_10kpc16384.dat"); # V1 time, V2 signal
s20 = read.table("inputs/s20-gw_10kpc16384.dat"); # V1 time, V2 signal

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
#true_data = read.table("../simulation_data/s20_data_g2.dat", sep = ",", comment.char = "#",header = TRUE); # true values
true_data = read.table("inputs/s20_data_g2.dat", sep = ",", comment.char = "#",header = TRUE); # true values
colnames(true_data) = c ("time","x");

##################################
### Frequentist analysis model ###
##################################
# data to generate model
fits_data = read.table("../simulation_data/A-A_fits_data_g2.dat", sep = ",", comment.char = "#",header = TRUE);
fits_data = read.table("inputs/A-A_fits_data_g2.dat", sep = ",", comment.char = "#",header = TRUE);

colnames(fits_data) = c("r", "f");

#mod = lm(r~ 0 + f + I(f^2) + I(f^3), data = fits_data); # first model
#summary(mod);
mod = lm(r~ 0 + f + I(f^3), data = fits_data);# final model
summary(mod);


###########################################
### generate noisy data and add signal ####
###########################################
d = data_generator(fs, duration, wvf.df, ampl, fcut, actPlot=FALSE)
noisydata = data.frame("V1"=d$t,"V2"=d$y)

# plot a spectrogram
r = specPdgrm(d$y,d$t, l=200, p=90, fs=fs, actPlot = TRUE, 
              logPow = T, zoomFreq=c(0,1)); # generating the spectrogram

############################
### FREQUENTIST ANALYSIS ###
############################

# find gmode frequency and compute the overlap with the model
covpbb(noisydata, mod=mod, l=200, p=90, fs=fs,um=10,dm=10, 
       thruth_data=true_data, actPlot=TRUE)


covpbb1(noisydata, mod=mod, l=200, p=90, fs=fs,um=10,dm=10, 
        thruth_data=true_data, actPlot=TRUE,
        limFreq = c(1000, 1200, 1300, 1400, Inf));

R = repcovpbb(wvf.df, duration, ampl = c(10), fcut,
              mod, N = 10, snr = NULL, movGmode = 11, um = 3, dm = 3,
              movBand = 5, timeGmode = NULL, l=200, p=90, fs=fs, limFreq=1000,
              thruth_data=true_data)

boxplot(R$covpbb,ylim = c(0, 1),xlab='amplitude factor',ylab='ratio',main='Frequentist')

#########################
### BAYESIAN ANALYSIS ###
#########################

postSamples = read.table("postsamples.txt", header = T); # importing posterios samples

covpbbBayes(noisydata, postSamples, l=200, p=90, fs=fs, movGmode = 11, 
            um=10, dm=10, movBand=5, timeGmode = NULL, 
            thruth_data=true_data, actPlot =TRUE, fmean = mean(fits_data$f))

#########################
### Multiple Analysis ###
#########################

# Multiple Analysis for different distances and frequency threshold

N    = 100;          # number of replications
dist = c(1,2,3,4,5); # distances = 10/amplitudFactor
out1 = out2 = NULL;  # to store outputs

for(d in dist){
  print(d)
  for(i in 1:N){
    #print(i)
    set.seed(i);
    newdata = data_generator(fs, duration, wvf.df, ampl = 10/d, fcut, actPlot=FALSE);
    noisydata = data.frame("V1"=newdata$t, "V2"=newdata$y);
    
    x = covpbb1(noisydata, mod=mod, l=200, p=90, fs=fs,um=10,dm=10, 
                thruth_data=true_data, actPlot=FALSE,
                limFreq = c(1000, 1100, 1200, Inf));
    
    out1 = rbind(out1, cbind(d, x$covpbb));   # coverage probability 
    out2 = rbind(out2, cbind(d, x$residual)); # statistics of residuals  
  }
}

# Boxplots # Note that there are NA values in out1
boxplot(covpbb~d+limFreq, out1, xlab = "distance : freq threshold");
boxplot(covpbb~limFreq+d, out1, xlab = "freq threshold : distance");
