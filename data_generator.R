library ("stats")
library ("signal")
library ("seewave")
library ("psd")
library ("pracma")


########################################################################
signal_generator = function(name, fs, actPlot=TRUE, pbOff=TRUE){
########################################################################
  
  # name: s11.2--LS220--GravA
  folder="inputs/New_simulations/"
  
  if (name=="s20-gw_10kpc16384"){
    gw_filename = paste(folder, name, ".dat", sep="")
    truedata_filename=paste(folder, "Ratios/s20_data_g2.dat", sep="")
  }else{
    gw_filename = paste(folder, "gw_", name, ".dat", sep="")
    truedata_filename = paste(folder, "Ratios/", name, "_mass_r.csv", sep="")
  }
  
  metadata_filename = paste(folder,"metadata.csv", sep="")
  
  if (actPlot==TRUE){
    print(gw_filename)
    print(truedata_filename)
    print(metadata_filename)
  }
  # Signal
  sXX = read.table(gw_filename); # V1 time, V2 signal
  colnames(sXX) = c ("time","hoft");
  duration=length(sXX$time)/16384
  
  # True data to define the ratio
  # For g2 modes, ratio (x variable in TF19) is x = Mpns / Rpns^2 
  true_data = read.table(truedata_filename, sep = ",", comment.char = "#",header = TRUE);
  
  if (name!="s20-gw_10kpc16384"){
    true_data = cbind(true_data$time, true_data$mass_pns / true_data$r_pns^2);
  }
  colnames(true_data) = c ("time", "ratio");
  true_data = as.data.frame(true_data);
  
  # Metadata
  meda = read.csv(metadata_filename)
  colnames(meda) = c("name","tb")
  t_bounce=meda$tb[which(meda$name == name)]
  
  if (actPlot==TRUE){
    print(c("Number of samples at 16384 Hz:", length(sXX$time)))
    print(c("Waveform duration:", duration, "s"))
    print(paste(name,"tbounce=", t_bounce," s"))
  }
  
  # shift time such that t_bounce=0
  sXX$time=sXX$time-t_bounce  
  true_data$time=true_data$time-t_bounce 
  
  # signal sampled at 16384 Hz. Resampling at fs
  resamp_factor=16384/fs
  signaly=resample(sXX$hoft,1/16384,1/fs)
  
  # decimate time vector
  signalx=seq(1,(fs*duration),by=1)
  for (i in 1:(fs*duration)) {
    signalx[i]=mean(sXX$time[((i-1)*resamp_factor+1):((i-1)*resamp_factor+resamp_factor)])
  }
  sYY = data.frame("time"=signalx,"hoft"=signaly)
  
  # remove times corresponding to the post-bounce period that is removed (100 ms) for all wvfs
  if (pbOff==TRUE){
    t_start=t_bounce+0.100
    true_data=subset(true_data, true_data$time >= 0.100)
    sYY=subset(sYY, sYY$time >= 0.100)  
  }
  
  if (actPlot == TRUE){
    plot(sXX$time, sXX$hoft, type="l", xlab="Time after bounce [s]", ylab="h(t)",, xlim=c(0,duration),pch=1)
    points(sYY$time, sYY$hoft, col="red",pch=2)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    title(name)
    
    leg = c("wvf resampled", "first 100ms after bounce removed")
    col = c("black","red")
    legend (x=duration*.1,y=max(sXX$hoft)*.9,legend=leg,col=col,pch=c(1,2))
  }
  
  return (list(wvf=sYY, true_data=true_data, duration=duration))
  
}


########################################################################
data_generator = function (fs, duration, wvf.df, ampl=1, filtering, fcut, actPlot=TRUE){
########################################################################
  # fs: sampling frequency
  # duration: duration (in second) of the output time serie
  # wvf: dataframe (time=time, hoft=h(t)) that contains the signal waveform sampled at fs
  # ampl: multiplication factor to change the source distance
  # fcut: lower frequency cutoff of the high pass filter. If NULL no HP filter.
  # filtering: "HP" --> fcut must be set
  #            "whitening" --> fcut is not used
  #            default is "HP"
  # The wvf is centered in the noise vector
  # 
  # output: d$t: time vector
  #         d$x: noise+signal
  #         d$y: filtered (noise+signal)
  
  
  # check that duration is an integer > 0
  if ((duration%%1)>0){
    print("duration must be a positive integer")
    return()
  }
 
  n=duration*fs
  
 
  n_0=n
  wvf_opt=0
  print(length(wvf.df))
  # Prepration in case a signal is added  
  if (length(wvf.df)!=0) {
    wvf_opt=1
    wvf_size=length(wvf.df$hoft)
    if (wvf_size>n){
      print("The waveform length is larger than the requested data vector size. Please check")
      return(NULL)
    }

#    if (2*wvf_size>n){
#      print("The internal vector's duration is doubled to avoid filtering border effect.")
#      print("The output vector size remains.")
      n=2*duration*fs
#    }

    #if (wvf_size%%2>0){
    #  print("The waveform length should be an even number. Please check")
    #  return(NULL)
    #}
  }
  
  # Noise generation
  freq2 = fs*fftfreq(n)          # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:int(n/2)]        # one-sided frequency vector
  psd=aLIGO_PSD_new(freq2,2)     # two-sided PSD          

  T = seq(1, n, by = 1)  

  X = rnorm(n, mean=0, sd=1);           # Gaussian white noise
  XX = fft(X); # FFT computing and normalization
  XXX = XX*sqrt(psd);                   # Coloring
  Y = fft(XXX, inverse = TRUE);         # FFT inverse
  Y = Re(Y)/sqrt(fs);                   # noise in time domain

  X1 = rnorm(n, mean=0, sd=1);            # Gaussian white noise
  XX1 = fft(X1); # FFT computing and normalization
  XXX1 = XX1*sqrt(psd);                   # Coloring
  Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
  Y1 = Re(Y1)/sqrt(fs);                   # noise in time domain
  
  # Signal addition (centered at the middle of the data vector to avoid filtering leakage
  # at the beggining and end).
  if (wvf_opt==1){
    ind1=floor((n-wvf_size)/2.)
    
    for (i in 1:wvf_size){
      Y[ind1+i]=Y[ind1+i]+ampl*wvf.df$hoft[i]
    }
  }

  if (filtering == "HP"){
    if (is.null(fcut)) {
      stop("fcut is not defined")
    }
    # filter noise + signal 
    # filtfilt : zero phase filter (forward& backward)
    myfilter=butter(n=4,W=fcut/(fs/2),type="high")
    YY=filtfilt(filt=myfilter,x=Y)          
  }
  else 
    if (filtering == "whitening"){
      # whiten the data
      ar_model <- stats::ar(Y1,order.max=3, aic=FALSE ,method=c("yule-walker"), demean=TRUE);
      b <- stats::filter(x=Y, filt=c(1, -ar_model$ar[1], -ar_model$ar[2], -ar_model$ar[3]), method="convolution",  sides = 1);
      YY=b
    }
    else if (filtering == "whitening_bis"){
      # compute the PSD
      psdest <- pspectrum(Y1,Y1.frqsamp=fs, ntap.init=6, Nyquist.normalize = TRUE, plot=FALSE,verbose = FALSE)
      psdwhitening=sqrt(fs)*psdest$spec/(n*sqrt(n))
      
      a = sqrt(fs) * fft(Y) / (n*sqrt(n)); # FFT computing and normalization
      b = a/sqrt(psdwhitening)             # whitening
      c = fft(b, inverse = TRUE);           # FFT inverse
      YY = Re(c)*sqrt(n);
      myfilter=butter(n=4,W=10/(fs/2),type="high")
      YY=filtfilt(filt=myfilter,x=YY)          
    }
  else {
    print("No filtering method specify")
    YY=Y
  }

  # select the orginal data size
  Tf = seq(1, n_0, by = 1)

  for (i in 1:n_0){
    Tf[i]=wvf.df$time[1]+i/fs
  }
  
  Yf = seq(1, n_0, by = 1)
  YYf = seq(1, n_0, by = 1)

  for (i in 1:n_0){
    Yf[i]=Y[ind1+i]
    YYf[i]=YY[ind1+i]
  }

  if (actPlot==TRUE){
    if (filtering == "HP"){
      ts.plot(Y); # noise only
      points(T, Y, col="black", type="l", pch=1, panel.first = grid())
      points(T, YY, col="red", type="l", pch=2);                                # (noise + signal) filtered
      T_wvf=seq(ind1,ind1+wvf_size-1,by=1)
      points(T_wvf,(wvf.df$hoft)*ampl,col="green",type="l",pch=3);  # signal only

      leg = c("noise", "(noise+signal) HP filtered", "signal only")
      col = c("black","red","green")
      legend (x=0,y=max(Y)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))
    }

    # Fourier transform
    YYFT = sqrt(fs) * fft(YY) / (n*sqrt(n)); # FFT computing and normalization
    plot (freq1, abs(YYFT)[1:int(n/2)], log="xy", type="l", xlab="frequency",ylab="ASD", col="grey",xlim=c(1, 4*10^3),panel.first = grid())
    legend (x=10, y=min(abs(YYFT))*1.2, legend="signal+noise FT")
 
    # Output vectors time series

    if (filtering == "HP"){
      plot (Tf, Yf, type="l", col="black")
      points (Tf, YYf, type="l", col="red")
      legend (x=Tf[1], y=max(Yf)*.9, col=c("black","red"), legend=c("noise+signal", "noise+signal after HP filter"))
    }else{
      plot (Tf, Yf, type="l", col="black")
      legend (x=Tf[1], y=max(Yf)*.9, legend="noise+signal")
    }
  }

  return(list(t=Tf,x=Yf,y=YYf))
  }

########################################################################
data_generator1 = function (fs, duration, wvf.df, ampl=1, filtering, setseed=0, actPlot=TRUE, verbose=FALSE){
########################################################################
  # fs: sampling frequency
  # duration: duration (in second) of the output time serie
  # wvf: dataframe (time=time, hoft=h(t)) that contains the signal waveform sampled at fs
  # ampl: multiplication factor to change the source distance
  # fcut: lower frequency cutoff of the high pass filter. If NULL no HP filter.
  # filtering method:
  #   "HP" : The fcut parameter is fixed internally (15 Hz)
  #   "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #   "AR" : AR model
  #   "prewhiten": use the R prewhiten function
  # 
  # setseed: if null, random seed. Otherwise set to the value
  # The wvf is centered in the noise vector
  # 
  # output: d$t: time vector
  #         d$x: noise+signal
  #         d$y: filtered (noise+signal)
  
  
  # check that duration is an integer > 0

  if ((duration%%1)>0){
    print("duration must be a positive integer")
    return()
  }
 
  n=duration*fs 
  if (verbose==TRUE){
    print(c("Size of the vectors: ", n))
    print(length(wvf.df))
  }
  
  n_0=n
  wvf_opt=0

  # Prepration in case a signal is added  
  if (length(wvf.df)!=0) {
    wvf_opt=1
    wvf_size=length(wvf.df$hoft)
    if (wvf_size>n){
      print("The waveform length is larger than the requested data vector size. Please check")
      return(NULL)
    }

    #if (wvf_size%%2>0){
    #  print("The waveform length should be an even number. Please check")
    #  return(NULL)
    #}

    # double the vector size
    duration=2*duration
    n=duration*fs
  }
  
  # Noise generation
  freq2 = fs*fftfreq(n)          # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:int(n/2)]        # one-sided frequency vector
  psd=aLIGO_PSD_new(freq2,2)     # two-sided PSD          

  # compute noise sigma
  s0 <- sqrt(2*trapz(freq1,psd[1:int(n/2)]))
  if (verbose==TRUE){print(c("aLIGO ASD noise rms:", s0))}
  
  T = seq(1, n, by = 1)  

  if (setseed>0){set.seed(setseed)}
  
  X = rnorm(n, mean=0, sd=1);           # Gaussian white noise
  XX = fft(X);                          # FFT computing and normalization
  XXX = XX*sqrt(psd);                   # Coloring
  Y = fft(XXX, inverse = TRUE);         # FFT inverse
  Y = Re(Y)/sqrt(fs);                   # noise in time domain

  # Signal addition (centered at the middle of the data vector to avoid filtering leakage
  # at the beggining and end).
  if (wvf_opt==1){
    ind1=floor((n-wvf_size)/2)
    
    for (i in 1:wvf_size){
      Y[ind1+i]=Y[ind1+i]+ampl*wvf.df$hoft[i]
    }
  }

  if (filtering == "HP"){
    fcut=15
    # filtfilt : zero phase filter (forward& backward)
    myfilter=butter(n=4,W=fcut/(fs/2),type="high")
    YY=filtfilt(filt=myfilter,x=Y)          
  }
  else 
    if (filtering == "AR"){
      # generate another noise TS
      X1 = rnorm(n, mean=0, sd=1);            # Gaussian white noise
      XX1 = fft(X1);                          # FFT computing
      XXX1 = XX1*sqrt(psd);                   # Coloring
      Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
      Y1 = Re(Y1)/sqrt(fs);                   # noise in time domain

      ar_model <- stats::ar(Y1,order.max=10, aic=FALSE ,method=c("yule-walker"), demean=TRUE);
      b <- stats::filter(x=Y, filt=c(1, -ar_model$ar[1], -ar_model$ar[2], -ar_model$ar[3], -ar_model$ar[4], -ar_model$ar[5], -ar_model$ar[6], -ar_model$ar[7], -ar_model$ar[8], -ar_model$ar[9], -ar_model$ar[10]), method="convolution",  sides = 1);
      b[1]=b[2]=b[3]=b[4]=b[5]=b[6]=b[7]=b[8]=b[9]=b[10]=b[11]
      YY=b
    }
    else
      if (filtering == "spectrum"){
        # generate another noise TS
        X1 = rnorm(n, mean=0, sd=1);           # Gaussian white noise
        XX1 = fft(X1);                          # FFT computing
        XXX1 = XX1*sqrt(psd);                   # Coloring
        Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
        Y1 = Re(Y1)/sqrt(fs);                   # noise in time domain

        # compute the PSD
        psdest <- pspectrum(Y1,Y1.frqsamp=fs, ntap.init=6, Nyquist.normalize = TRUE, plot=FALSE,verbose = FALSE)
        psdwhitening=rep(0, n);
        for(i in 1:(int(n/2))){
          psdwhitening[i]=psdest$spec[i]
          psdwhitening[n+1-i]=psdest$spec[i]
        }

        psdwhitening=psdwhitening*fs/2./n/n
    
        a = fft(Y)                        # FFT computing and normalization
        b = a/sqrt(psdwhitening)          # whitening
        c = fft(b, inverse = TRUE);       # FFT inverse
        YY = s0/10*Re(c)/n/sqrt(fs);         # Normalisation factor of the 2 FFTs

      }
      else 
        if (filtering == "prewhiten"){
          # prewhiten
          myts <- ts(Y, start=0, end=duration, frequency=fs)
          #myts <- ts(Y)
          myts <- prewhiten(myts, AR.max=100, zero.pad="rear", plot=FALSE, verbose=FALSE)

          YY <- myts[['prew_ar']][1:n]
        }
        else {
          print("No filtering method specify")
          YY=Y
        }


  # select the original data size
  Tf = seq(1, n_0, by = 1)

  for (i in 1:n_0){
    Tf[i]=wvf.df$time[1]+i/fs
  }
  
  Yf = seq(1, n_0, by = 1)
  YYf = seq(1, n_0, by = 1)

  for (i in 1:n_0){
    Yf[i]=Y[ind1+i]
    YYf[i]=YY[ind1+i]
  }

  if (actPlot==TRUE){
    if (filtering == "HP" || filtering == "spectrum" || filtering == "prewhiten" || filtering == "AR"){
      ts.plot(Y); # noise only
      points(T, Y, col="black", type="l", pch=1, panel.first = grid())
      points(T, YY, col="red", type="l", pch=2);                                # (noise + signal) filtered
      T_wvf=seq(ind1,ind1+wvf_size-1,by=1)
      points(T_wvf,(wvf.df$hoft)*ampl,col="green",type="l",pch=3);  # signal only

      leg = c("noise", "(noise+signal) filtered", "signal only")
      col = c("black","red","green")
      legend (x=0,y=max(Y)*.9,legend=leg,cex=.8,col=col,pch=c(1,3))
    }

    # Fourier transform
    YYFT = sqrt(fs) * fft(YY) / (n*sqrt(n)); # FFT computing and normalization
    plot (freq1, abs(YYFT)[1:int(n/2)], log="xy", type="l", xlab="frequency",ylab="ASD", col="grey",xlim=c(1, 4*10^3),panel.first = grid())
    legend (x=10, y=min(abs(YYFT))*1.2, legend="signal+noise FT")
 
    # Output vectors time series

    if (filtering == "HP" || filtering == "spectrum" || filtering == "prewhiten" || filtering == "AR"){
      plot (Tf, Yf, type="l", col="black")
      points (Tf, YYf, type="l", col="red")
      legend (x=Tf[1], y=max(Yf)*.9, col=c("black","red"), legend=c("noise+signal", "noise+signal after filter"))
    }else{
      plot (Tf, Yf, type="l", col="black")
      legend (x=Tf[1], y=max(Yf)*.9, legend="noise+signal")
    }
  }

  return(list(t=Tf,x=Yf,y=YYf))
  }
  
  
########################################################################
noise_generator = function (fs, duration, filtering, setseed=0, actPlot=TRUE, verbose=FALSE){
########################################################################
  # fs: sampling frequency
  # duration: duration (in second) of the output time serie
  # filtering method:
  #   "HP" : The fcut parameter is fixed internally (15 Hz)
  #   "spectrum" : the data are whiten in Fourier domain using the noise spectrum estimate
  #   "AR" : AR model
  #   "prewhiten": use the R prewhiten function
  # setseed: if null, random seed. Otherwise set to the value
  #
  # output: d$t: time vector
  #         d$x: noise
  
  # check that duration is an integer > 0
  
  if ((duration%%1)>0){
    print("duration must be a positive integer")
    return()
  }
  
  n=20*duration*fs 
  if (verbose==TRUE){print(c("Size of the vectors: ", n))}
  
  # Noise generation
  freq2 = fs*fftfreq(n)          # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:int(n/2)]        # one-sided frequency vector
  psd=aLIGO_PSD_new(freq2,2)     # two-sided PSD          
  
  #psd=PSD_3G(freq2,2,'ET_B',verbose)

  # compute noise sigma
  s0 <- sqrt(2*trapz(freq1,psd[1:int(n/2)]))
  if (verbose==TRUE){print(c("aLIGO ASD noise rms:", s0))}
    
  T = seq(1, n, by = 1)
  
  if (setseed>0){set.seed(setseed)}

  X = rnorm(n, mean=0, sd=1);           # Gaussian white noise
  XX = fft(X);                          # FFT computing
  XXX = XX*sqrt(psd);                   # Coloring
  Y = fft(XXX, inverse = TRUE);         # FFT inverse
  Y = Re(Y)/sqrt(fs);                   # noise in time domain
  
  if (filtering == "HP"){
    fcut=15
    # filtfilt : zero phase filter (forward& backward)
    myfilter=butter(n=4,W=fcut/(fs/2),type="high")
    YY=filtfilt(filt=myfilter,x=Y)          
  }
  else 
    if (filtering == "AR"){
      # generate another noise TS
      X1 = rnorm(n, mean=0, sd=1);            # Gaussian white noise
      XX1 = fft(X1);                          # FFT computing
      XXX1 = XX1*sqrt(psd);                   # Coloring
      Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
      Y1 = Re(Y1)/sqrt(fs);                   # noise in time domain

      ar_model <- stats::ar(Y1,order.max=10, aic=FALSE ,method=c("yule-walker"), demean=TRUE);
      b <- stats::filter(x=Y, filt=c(1, -ar_model$ar[1], -ar_model$ar[2], -ar_model$ar[3], -ar_model$ar[4], -ar_model$ar[5], -ar_model$ar[6], -ar_model$ar[7], -ar_model$ar[8], -ar_model$ar[9], -ar_model$ar[10]), method="convolution",  sides = 1);
      b[1]=b[2]=b[3]=b[4]=b[5]=b[6]=b[7]=b[8]=b[9]=b[10]=b[11]
      YY=b
    }
    else 
      if (filtering == "spectrum"){
        # generate another noise TS
        X1 = rnorm(n, mean=0, sd=1);           # Gaussian white noise
        XX1 = fft(X1);                          # FFT computing
        XXX1 = XX1*sqrt(psd);                   # Coloring
        Y1 = fft(XXX1, inverse = TRUE);         # FFT inverse
        Y1 = Re(Y1)/sqrt(fs);                   # noise in time domain

        # compute the PSD
        # 
        #myts <- ts(Y1, start=0, end=duration, frequency=fs)
        psdest <- pspectrum(Y1,Y1.frqsamp=fs, ntap.init=6, Nyquist.normalize = TRUE, plot=FALSE,verbose = FALSE)
        psdwhitening=rep(0, n);
        for(i in 1:(int(n/2))){
          psdwhitening[i]=psdest$spec[i]
          psdwhitening[n+1-i]=psdest$spec[i]
        }

        psdwhitening=psdwhitening*fs/2./n/n
    
        a = fft(Y)                        # FFT computing and normalization
        b = a/sqrt(psdwhitening)          # whitening
        c = fft(b, inverse = TRUE);       # FFT inverse
        YY = s0/10*Re(c)/n/sqrt(fs);      # Normalisation factor of the 2 FFTs

        #    myfilter=butter(n=4,W=10/(fs/2),type="high")
        #    YY=filtfilt(filt=myfilter,x=YY)          
      }
      else 
        if (filtering == "prewhiten"){
          #fcut=10
          #filtfilt : zero phase filter (forward& backward)
          #myfilter=butter(n=4,W=fcut/(fs/2),type="high")
          #Y=filtfilt(filt=myfilter,x=Y)     

          # prewhiten
          myts <- ts(Y, start=0, end=duration, frequency=fs)
          #myts <- ts(Y)
          myts <- prewhiten(myts, AR.max=100, zero.pad="rear", plot=FALSE, verbose=FALSE)
          
          YY <- myts[['prew_ar']][1:n]
          
        }
        else {
          print("No filtering method specify")
          YY=Y
        }
  
  # select the orginal data size
  Tf = seq(1, n, by = 1)
  
  for (i in 1:n){
    Tf[i]=i/fs
  }
  

  if (actPlot==TRUE){
    ts.plot(Y); # noise only
    points(T, Y, col="black", type="l", pch=1, panel.first = grid())
    points(T, YY, col="red", type="l",pch=2)
      
    leg = c("raw noise","prewhiten noise")
    col = c("black","red")
    legend (x=0,y=max(Y)*.9,legend=leg,cex=.8,col=col,pch=c(1,2))
    
    # spectrum estimated
    psdest <- pspectrum(Y, Y.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, plot=FALSE,verbose = FALSE)
    psdest_filtered <- pspectrum(YY, YY.frqsamp=fs, ntap.init=NULL, Nyquist.normalize = TRUE, plot=FALSE,verbose = FALSE)
 
    # Fourier transform
    YFT = sqrt(fs) * fft(Y) / (n*sqrt(n));
    WFT = sqrt(fs) * fft(YY) / (n*sqrt(n));
    plot (freq1, abs(YFT)[1:int(n/2)], log="xy", type="l", xlab="frequency", ylab="ASD", col="grey", xlim=c(1, fs/2), ylim=c(1e-26,1e-20), pch=1, panel.first = grid())

    lines(freq1, sqrt(psd[1:int(n/2)]), col="red", pch=3)
    lines(fs*psdest$freq, sqrt(psdest$spec)*sqrt(fs/2)/n, col="blue", pch=2)
   
    lines(freq1, abs(WFT)[1:int(n/2)], col="black", pch=4)
    lines(fs*psdest_filtered$freq[1:int(n/2)], sqrt(psdest_filtered$spec[1:int(n/2)])*sqrt(fs/2)/n, col="green", pch=5)

    legend_str=c("col noise FT", "col noise spectrun", "ASD model", "filtered FT", "filtered spectrum")
    legend (x=100, y=min(abs(tail(YFT,-1)))*50000, legend=legend_str, col=c("grey","blue","red","black","green"), pch=c(1,2,3,4,5))   

    s1 <- sqrt(2*trapz(fs*psdest$freq[1:int(n/2)], psdest$spec[1:int(n/2)]*fs/(2*n*n)))
    if (verbose==TRUE){print(c("colored noise rms:", s1))}
    
    s2 <- sqrt(2*trapz(fs*psdest_filtered$freq[1:int(n/2)], psdest_filtered$spec[1:int(n/2)]*fs/(2*n*n)))
    if (verbose==TRUE){print(c("filtered noise rms:", s2))}
    
    # Output vectors time series
    
#    if (filtering == "HP"){
#      plot (Tf, Yf, type="l", col="black")
#      points (Tf, YYf, type="l", col="red")
#      legend (x=Tf[1], y=max(Yf)*.9, col=c("black","red"), legend=c("noise+signal", "noise+signal after HP filter"))
#    }
  }
  
  return(list(t=Tf,x=Y))
}


########################################################################
fftfreq = function(n, d = 1){
########################################################################
  # surogate for the numpy fft.fftfreq function that generates the two sided 
  # frequency vector. Defaults d=1 means sampling frequency is 1. 
  # https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
  #
  # n: samples number
  # d: sample spacing (inverse of the sampling rate). Defaults to 1
  
  if(n%%2 == 0){# n is even
    out = c(seq(0, n/2-1, by = 1), seq(-n/2, -1, by=1)) / (d*n);
  }else{ # n is odd
    out = c(seq(0, (n-1)/2, by = 1), seq(-(n-1)/2, -1, by=1)) / (d*n);
  }
  
  return(out);
}

########################################################################
int = function(n){
########################################################################
  # https://stackoverflow.com/questions/31036098/what-is-the-difference-between-int-and-floor-in-python-3
  if(n < 0 ){
    return(-floor(abs(n)))
  }else{
    return(floor(n))
  }
}

########################################################################
aLIGO_PSD = function(f,type){
########################################################################
  # Original aLIGO PSD function used by Patricio
  cutoff = -109.35 + log(2e10);
  fn     = length(f);
  logpsd = rep(0, fn);
  if(f[1]==0){
    f[1]=f[2]
  }
  if(type == 1){
    for(i in 1:fn){
      x  = f[i]/215;
      x2 = x*x;
      logpsd[i] = log(1e-49) + log(x^(-4.14) -5/x2 + 111*(1-x2+0.5*x2*x2)/(1.+0.5*x2));
      
      if(logpsd[i]>cutoff){
        logpsd[i]=cutoff;
      }
    }
    output=exp(logpsd);
    
  }else{
    for(i in 1:(int(fn/2)+1)){
      x  = abs(f[i]/215);
      x2 = x*x;
      logpsd[i] = log(1e-49) + log(x^(-4.14) -5./x2 + 111.*(1-x2+0.5*x2*x2)/(1.+0.5*x2));
      
      if(logpsd[i]>cutoff){
        logpsd[i]=cutoff;
      }
      if(i>0){
        logpsd[fn-i]=logpsd[i];
      }
    }
    output=exp(logpsd)/2;            # Two sided PSD
  }
  return (output)
}

########################################################################
aLIGO_PSD_new = function(f,type){
########################################################################
  # aLIGO sensitivity curve: fit the data point from https://dcc.ligo.org/LIGO-T1800044/public
  # Type=1 --> one-sided PSD.
  # Type=2 --> two-sided PSD. 
  
  S1 = 5.0e-26;
  S2 = 1.0e-40;
  S3 = 1.4e-46;
  S4 = 2.7e-51;
  fcut = 10;
  cutoff=1e-42;
  fn     = length(f);
  output = rep(0, fn); #np.zeros(len(f))
  
  # to avoid issue with f=0
  if(f[1]==0){
    f[1]=f[2];
  }
  if(type == 1){
    for(i in 1:fn){
      x = abs(f[i]);
      output[i] =  S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }
    }
  }else{
    for(i in 1:(int(fn/2)+1)){ # range(int(len(f)/2)+1)
      x = abs(f[i]);
      output[i] = S1/(x^20) + S2/(x^4.05) + S3/(x^.5) + S4*((x/fcut)^2);
      if(output[i]>cutoff){
        output[i]=cutoff
      }

      # Wraparound frequency: f=0 is the first element (i=1), 
      # and all elements are symetric around index fn/2+1
      if(i>1){
        output[fn+2-i]=output[i]
      }
    }  
    output=output/2;            # Two sided PSD
#    output=shifter(output,-1)  # No more needed after correction 4 lines above
  }
  return(output);
}

shifter = function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}    

########################################################################
PSD_3G = function(f,type,detector,verbose=FALSE){
  ########################################################################
  # 3G sensitivity curves ref:
  # f: frequency vector
  # type=1 --> one-sided PSD.
  # type=2 --> two-sided PSD.
  # detector: name of the 3G detector
  
  psd_filename = sprintf("PSD/%s_data.txt", detector)
  #psd_filename = sprintf("PSD/%s_sensitivity.txt", detector)
  data = read.table(psd_filename);
  fs=4096
  cutoff=1e-42
  
  n = length(f)  
  fmin=f[1]
  if (type==1){
    fmax=f[n]
  } else{
    fmax=abs(f[n/2+1])}

  yl=data$V2[1]
  yr=data$V2[length(data$V2)]
  
  asd_func = approxfun(x = data$V1, y = data$V2, method = "linear",
                yleft=yl, yright=yr, rule = 1, f = 0, ties = "mean");
  
  if (type==1){
    asd = asd_func(f)
    psd = asd*asd
  }else{
    asd = rep(0, n);
    asd_1sided = asd_func(abs(f[1:int(n/2)]))
  
    if (length(asd_1sided) != int(n/2)){
      print ("Warning: ASD size is different from what is expected")
    }
    
    for(i in 1:(int(n/2))){ 
      asd[i]=asd_1sided[i];
      
      # Wraparound frequency: f=0 is the first element (i=1), 
      # and all elements are symetric around index n/2+1
      if(i>1){
        asd[n+2-i]=asd[i]
      }
    }  
    asd[n/2+1]=asd_func(abs(f[int(n/2)+1]))

    # Two sided psd
    asd=asd/2;
    psd=asd*asd;
  }
  
  for (i in 1:n){
    if (psd[i]>cutoff){
      psd[i]=cutoff
    }
  }
  
  if (verbose==TRUE){
    if (type==1){
      plot(f,psd,log="xy",col="blue",xlim=c(1, fs/2),pch=2)
      points(data$V1,data$V2*data$V2,col="red",type="l",pch=1)
    }else{
      plot(f,psd,log="y",col="blue",xlim=c(1, fs/2),pch=2)
      points(data$V1,0.25*data$V2*data$V2,col="red",type="l",pch=1)
    }
    leg = c(detector,"interpolated")
    col = c("red","blue")
    legend (x=500,y=psd[1]*0.8,legend=leg,cex=.8,col=col,pch=c(1,2))
  }
  

  return(psd)
}

########################################################################
compute_SNR = function(name, fcut){
########################################################################
  fs=4096
  signal=signal_generator(name, fs, actPlot=TRUE, pbOff=TRUE)

  n=length(signal$wvf$hoft)
  n2=2*fs*signal$duration        # zeropadding

  #print(c(n,n2))
  
  freq2 = fs*fftfreq(n2)         # two-sided frequency vector
  freq2[1]=0.001                 # to avoid plotting pb in logscale
  freq1=freq2[1:int(n2/2)]       # one-sided frequency vector
  psd=aLIGO_PSD_new(freq1,1)     # one-sided PSD          

  waveform=signal$wvf
  vec=rep(0,n2)
  for (i in 1:n){
    vec[n2/4+i]=vec[n2/4+i]+waveform$hoft[i]
  }  
  #w=hanning(n2)
  #vec=waveform$hoft
  
  hf=fft(vec);
  hf=hf[1:(n2/2)]
  
  hf=subset(hf,freq1-fcut<0)
  psd=subset(psd,freq1-fcut<0)
  freq1=subset(freq1, freq1-fcut<0)
  
  integrand=abs(hf*Conj(hf))/(fs*fs)
  snr=sqrt(4*trapz(freq1,integrand/psd))
  
  print(c(name,"SNR:",snr))
}
