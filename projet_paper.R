source("specPdgrm.R")
source("data_simulation.R")
source("functions.R")

library(signal)

#########################
### Statistical model ###
#########################
fits_data = read.table("inputs/A-A_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("inputs/CoCo_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("inputs/fits_data.dat", sep = ",");# data to generate model
colnames(fits_data) = c("r", "f");

# Linear model
mod = lm(r~ 0 + f + I(f^3), data = fits_data);# final model
summary(mod);

# Variable variance linear model
mu1 = "~ f - 1";
mu2 = "~ f + I(f^2) - 1";
mu3 = "~ f + I(f^2) + I(f^3) - 1";
mu4 = "~ f + I(f^3) - 1";
s1 = "~ f - 1"; 
s2 = "~ f + I(f^2) - 1";

Xm  = model.matrix(eval(parse(text=eval(parse(text=mu3)))), fits_data);
Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = TRUE);

############################
### Algorithm parameters ###
############################
# g-mode estimate: starting from the right, left or median

#gmode = c("left", "right", "median");
gmode = c("left");

#############################
### Simulation parameters ###
#############################
#detectors=c("aLIGO", "ALIGO", "CE1", "CE2", "ET_B", "ET_C", "ET_D")
signals=c("s20.0--SFHo", "s11.2--LS220", "s15.0--GShen", "s15.0--LS220", 
          "s15.0--SFHo", "s20.0--LS220", "s25.0--LS220", "s40.0--LS220")
detectors=c("aLIGO", "CE1", "CE2", "ET_B", "ET_C", "ET_D")
freq_min=c(200,200,200,200,200,200,200,200)
distance_step=c(.5,3,3,3,3,3)
signals=c("s20.0--SFHo")
detectors=c("aLIGO")
fs=4096
filtering_method="prewhiten"

# loop over N generation of noisy data and add signal

N=100
dist_nb=38
result<-matrix(nrow=N*dist_nb,ncol=6)

isig=0
for (signal_name in signals){
  isig=isig+1
  signal = signal_generator(fs, signal_name, actPlot=FALSE, verbose=FALSE)
  wvf=signal$wvf
  true_data=signal$true_data
  duration=signal$duration
  idet=0  
  for (detector in detectors){
    # To use always the same random noise for all waveforms % detectors
    set.seed(1)  
    idet=idet+1
    for (j in 1:dist_nb){
      dist = 1+(j-1)*distance_step[idet]
      for (i in 1:N){
        d=data_generator(fs, duration, wvf, ampl=10/dist, detector, 
                   filter = filtering_method, 
                   setseed=0, actPlot=FALSE, verbose=FALSE);
        noisydata = data.frame("V1"=d$t,"V2"=d$y);
    
        out = covpbb(noisydata, mod=fit, l=200, p=90, fs=fs,
                  um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(freq_min[isig], 500),
                  um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(1000, 1700),
                  gmode = gmode, thruth_data=true_data, actPlot=FALSE,
                  limFreq = c(1000));
    
        result[i+(j-1)*N,1]=dist
        result[i+(j-1)*N,2]=out$covpbb[1,1]
        result[i+(j-1)*N,3]=out$covpbb[1,2]
        result[i+(j-1)*N,4]=out$residual[1,1]
        result[i+(j-1)*N,5]=out$residual[1,2]
        result[i+(j-1)*N,6]=out$residual[1,3]
        
      }
      #print(out1$covpbb[2])
      ind1=1+(j-1)*N
      ind2=N+(j-1)*N
      print(sprintf("%s and signal %s @ distance: %f kpc. Covpbb mean:%f. Covpbb median: %f",
                    detector, signal_name, dist, mean(result[ind1:ind2,2]), median(result[ind1:ind2,2])))
    }
    filename=sprintf("results_AA_%s_f2_%s_%s.txt", filtering_method, signal_name, detector)
    write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)
  }
}
  
