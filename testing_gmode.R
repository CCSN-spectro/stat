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

# For g2 modes, ratio (x variable in TF19) is x = Mpns / Rpns^2 
true_data = read.table("inputs/s20_data_g2.dat", sep = ",", comment.char = "#",header = TRUE); # true values
colnames(true_data) = c ("time","x");

###############
### Testing ###
###############

# Defining time for the estimated g-modes
index = ints(n=dim(noisydata)[1], l=200, p=90); # from psplinePsd
# centred point for even "l"
mindx = index + cbind(rep(200/2 -1, dim(index)[1]), rep(-(200/2-1), dim(index)[1]));
# Ajusting data
timedata = seq(from=wvf.df$V1[1], to=wvf.df$V1[dim(noisydata)[1]], 
               length = dim(noisydata)[1]);
# mean time (centred point) for our g-mode estimates
timefreq = apply(mindx, 1, function(x) mean(timedata[c(x[1], x[2])]) );

###

i = 8;
d = 2;
#for(i in 51:80){
set.seed(i);
newdata = data_generator(fs, duration, wvf.df, ampl = 10/d, filtering="HP", fcut=fcut, actPlot=FALSE);
noisydata = data.frame("V1"=newdata$t, "V2"=newdata$y);
r = specPdgrm(noisydata$V2,noisydata$V1, l=200, p=90, fs=fs, actPlot = TRUE, 
              logPow = T, zoomFreq=c(0,1)); # generating the spectrogram
maxf = findGmodes(r, um_R = 0, dm_R = 8, m_R = 8, um_L = 8, dm_L = 0, m_L = 8,  
                  initfreq_R = c(800, 2000),initfreq_L = c(-Inf, Inf));

points(timefreq, maxf$maxf_med, col = "black");  # Median of the estimates 
points(timefreq, maxf$maxf_R, col = "blue"); # Starting from the right
points(timefreq, maxf$maxf_L, col = "red");  # Starting from the left 

#}

