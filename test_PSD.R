source("data_generator.R")

fs=4096
duration=1
fcut=15
n=duration*fs

freq2 = fs*fftfreq(n)   # two-sided frequency vector
freq1=freq2[1:int(n/2)]        # one-sided frequency vector
psd1=aLIGO_PSD_new(freq1,1)    # one-sided PSD

plot(freq1,sqrt(psd1),log="xy",type="l",xlab="frequency",ylab="ASD",pch=1,ylim=c(1e-24,1e-21),xlim=c(1, 4*10^3),panel.first = grid())

### generate noisy data and signal buried ###
d = data_generator(fs, duration, wvf=NULL, ampl=0, fcut, actPlot=FALSE)

FT = sqrt(fs) * fft(d$x) / (n*sqrt(n)); # FFT computing and normalization
points(freq1,2*abs(FT)[1:int(n/2)],type="l",col="red",pch=2)

FT = sqrt(fs) * fft(d$y) / (n*sqrt(n)); # FFT computing and normalization
points(freq1,2*abs(FT)[1:int(n/2)],type="l",col="green",pch=3)

points(freq1,sqrt(psd1),type="l")

leg.txt <- c("aLIGO ASD", "random noise","filtered")
col=c("black","red","green")
legend(x=250,y=10^-21,legend=leg.txt,cex=.8,col=col,pch=c(1,2,3))

