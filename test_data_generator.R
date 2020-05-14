source("data_generator.R")

fs=4096
duration=1
ampl=10
fcut=15
n=duration*fs

### signal ###
s20 = read.table("inputs/s20-gw_10kpc16384.dat"); # V1 time, V2 signal

signaly=resample(s20$V2,1/16384,1/fs)
signalx=decimate(s20$V1,16384/fs)

wvf.df = data.frame("V1"=signalx,"V2"=signaly)

### generate noisy data and signal buried ###
d = data_generator1(fs, duration, wvf.df, ampl, filtering="prewhiten", actPlot=TRUE, verbose=TRUE)
sprintf("noise and signal is whitened. Standard deviation: %g", sqrt(var(d$y)))

d1 = data_generator1(fs, duration, wvf.df, ampl, filtering="spectrum", actPlot=TRUE, verbose=TRUE)
sprintf("noise and signal is whitened. Standard deviation: %g", sqrt(var(d1$y)))

plot(d$t,d$y, col="black", type="l")
points(d1$t,d1$y, col="red", type="l")


freq1 = fs*fftfreq(n)          # two-sided frequency vector
freq = freq1[1:int(n/2)]
wbFT = sqrt(fs)*fft(d$y)/(n*sqrt(n));

plot(freq,2*abs(wbFT)[1:int(n/2)], log="y",type="l", xlab="frequency",ylab="ASD", col="cyan",ylim=c(min(abs(wbFT)),5*max(abs(wbFT))))

dp = data_generator(fs, duration, wvf.df, ampl, filtering="whitening_bis", fcut=fcut, actPlot=TRUE)
sprintf("noise and signal is high passed. Standard deviation: %g", sqrt(var(dp$y)))

