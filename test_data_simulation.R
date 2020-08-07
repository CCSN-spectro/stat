source("data_simulation.R")

detector="aLIGO"
fs=4096
duration=2
factor=1
dist=1
mc=noise_generator(factor,fs, duration, detector, setseed=1, filter="prewhiten",
                      actPlot=TRUE, verbose=TRUE)




out=signal_generator(fs, "s20-gw", 
                     actPlot=TRUE, verbose=TRUE)

d=data_generator(fs, duration, wvf, ampl=10/dist, detector, filter = "prewhiten", setseed=1, 
                 actPlot=TRUE, verbose=TRUE);

noisydata = data.frame("V1"=d$t,"V2"=d$y);





######################################

X=out$wvf$hoft
n=length(X)

XFT = fft(X)/sqrt(n);
psdest <- pspectrum(X, X.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, 
                    plot=FALSE,verbose=FALSE)

freq2=fs*fftfreq(n)          # two-sided frequency vector
freq2[1]=0.001                 # to avoid plotting pb in logscale
freq1=freq2[1:int(n/2)]        # one-sided frequency vector

plot (freq1, sqrt(2)*abs(XFT)[1:int(n/2)], log="xy", type="l", xlab="Frequency", ylab="ASD", 
      col="grey", xlim=c(1, fs/2), ylim=c(1e-24,1e-20), pch=1, panel.first = grid())

lines(fs*psdest$freq, sqrt(psdest$spec), col="black", pch=2)

n1=10*n
Y=rep(0,n1)
for (i in 1:n){
  Y[n1/4+i]=X[i]
}
YFT = fft(Y)/sqrt(n);
psdest <- pspectrum(Y, Y.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, 
                    plot=FALSE,verbose=FALSE)

freq2=fs*fftfreq(n1)          # two-sided frequency vector
freq2[1]=0.001                 # to avoid plotting pb in logscale
freq1=freq2[1:int(n1/2)] 

points (freq1, sqrt(2)*abs(YFT)[1:int(n1/2)], log="xy", type="l", xlab="Frequency", ylab="ASD", 
      col="red", xlim=c(1, fs/2), ylim=c(1e-24,1e-20), pch=1, panel.first = grid())

lines(fs*psdest$freq, sqrt(psdest$spec), col="blue", pch=2)


detector="aLIGO"
fs=10000
duration=10

compute_SNR("s20-gw", 10)

n=fs*duration

freq2=fs*fftfreq(n)          # two-sided frequency vector
freq2[1]=0.001                 # to avoid plotting pb in logscale
freq1=freq2[1:int(n/2)]        # one-sided frequency vector

psd=aLIGO_PSD_new(freq2, 2)
set.seed(1)

X = rnorm(n, mean=0, sd=1);           # Gaussian white noise
XX = fft(X);                          # FFT computing
XXX = XX*sqrt(psd);                   # Coloring
Y = fft(XXX, inverse = TRUE);         # FFT inverse
Y = Re(Y)/n;                   # noise in time domain
print(Y[1:100])

T=seq(1,n,by=1)
plot(T,Y,col="black",type="l",panel.first=grid())

psdest <- pspectrum(Y, Y.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, 
                    plot=FALSE,verbose=FALSE)

YFT = fft(Y)/sqrt(n);

plot (freq1, 2*abs(YFT)[1:int(n/2)], log="xy", type="l", xlab="Frequency", ylab="ASD", 
      col="grey", xlim=c(1, fs/2), ylim=c(1e-24,1e-20), pch=1, panel.first = grid())

lines(fs*psdest$freq, sqrt(psdest$spec), col="blue", pch=2)
lines(freq1, sqrt(2*psd[1:int(n/2)]), col="red", pch=3) 






folder = "inputs/New_simulations/"
metadata_filename = paste(folder,"metadata.csv", sep="")
meda = read.csv(metadata_filename)
colnames(meda) = c("name","wvf_name", "ratio_name", "tb")
t_bounce=meda$tb[which(meda$name == name)]

# test of signal 
name="s20-gw_10kpc16384"
name="s11.2--LS220--GravA"
name="s15.0--GShen--GravA"
name="s15.0--LS220--GravA"
name="s15.0--SFHo--GravA"
name="s20.0--LS220--GravA"
name="s25.0--LS220--GravA"
name="s40.0--LS220"
fs=4096
signal = signal_generator(name, fs)


# test of noise

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

############ test of noise_generator
fs=4096
duration=1




