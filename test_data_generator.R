source("data_generator.R")

fs=4096
duration=1
ampl=5
fcut=15

### signal ###
s20 = read.table("inputs/s20-gw_10kpc16384.dat"); # V1 time, V2 signal

signaly=resample(s20$V2,1/16384,1/fs)
signalx=decimate(s20$V1,16384/fs)

wvf.df = data.frame("V1"=signalx,"V2"=signaly)

### generate noisy data and signal buried ###
d = data_generator(fs, duration, wvf.df, ampl, fcut, actPlot=TRUE)

