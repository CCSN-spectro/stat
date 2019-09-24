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

Obergaunlinguer_data = read.csv("inputs/Obergaunlinguer_dataNEW.csv"); # true values
true_ratios = Obergaunlinguer_data$Mpns / (Obergaunlinguer_data$Rpns^(2))

############################
### FREQUENTIST ANALYSIS ###
############################

#############
### Model ###
#############
fits_data = read.table("inputs/fits_data.dat", sep = ",");# data to generate model
colnames(fits_data) = c("f", "r");

#mod = lm(r~ 0 + f + I(f^2) + I(f^3), data = fits_data); # first model
#summary(mod);
mod = lm(r~ 0 + f + I(f^3), data = fits_data);# final model
summary(mod);


# generate noisy data and add signal
d = data_generator(fs, duration, wvf.df, ampl, fcut, actPlot=FALSE)
noisydata = data.frame("V1"=d$t,"V2"=d$y)

# plot a spectrogram
r = specPdgrm(d$y,d$t, l=200, p=90, fs=fs, actPlot = TRUE, 
              logPow = T, zoomFreq=c(0,1)); # generating the spectrogram

# find gmode frequency and compute the overlap with the model
covpbb(noisydata, mod=mod, l=200, p=90, fs=fs,um=10,dm=10, 
       Obergaunlinguer_data=Obergaunlinguer_data, actPlot=TRUE)





