source("specPdgrm.R")
source("data_generator.R")
source("functions.R")

library(signal)

##################
### parameters ###
##################
#name="s20-gw_10kpc16384"
#name="s11.2--LS220--GravA"
#name="s15.0--GShen--GravA"
name="s15.0--LS220--GravA"
#name="s15.0--SFHo--GravA"
#name="s20.0--LS220--GravA"
#name="s25.0--LS220--GravA"
#   name="s40.0--LS220"

fs=4096
ampl=10
filtering_method="prewhiten"

#############
### Model ###
#############
fits_data = read.table("inputs/A-A_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("inputs/CoCo_fits_data_g2.dat", sep = ",");# data to generate model
#fits_data = read.table("inputs/fits_data.dat", sep = ",");# data to generate model
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

Xm  = model.matrix(eval(parse(text=eval(parse(text=mu3)))), fits_data);
Xs  = model.matrix(eval(parse(text=eval(parse(text=s2)))), fits_data);

# Variable variance linear model
fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = TRUE);

#############################
### Signal and true ratio ###
#############################
signal = signal_generator(name, fs)

wvf.df=signal$wvf
true_data=signal$true_data
duration=signal$duration



############################
### FREQUENTIST ANALYSIS ###
############################

# g-mode estimate: starting from the right, left or median
#gmode = c("left", "right", "median");
gmode = c("left");
#gmode = c("right"); 

###########
### Run ###
###########

# loop over N generation of noisy data and add signal

actplot=TRUE
N=1
dist_nb=25
result<-matrix(nrow=N*dist_nb,ncol=12)

set.seed(1)

for (j in 1:dist_nb){
  dist = 1+(j-1)*.25
  print(c("distance:",dist,"kpc"))
  for (i in 1:N){
    d = data_generator1(fs, duration, wvf.df, ampl=10/dist, filtering = filtering_method, actPlot=FALSE);
    noisydata = data.frame("V1"=d$t,"V2"=d$y);
    out = covpbb(noisydata, mod=mod, l=200, p=90, fs=fs,
              um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(100, 500),
              um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(1000, 1700),
              gmode = gmode,
              thruth_data=true_data, actPlot=TRUE,
              limFreq = c(1000));

    out1 = covpbb(noisydata, mod=fit, l=200, p=90, fs=fs,
               um_L = 8, dm_L = 0, m_L = 8, initfreq_L = c(100, 500),
               um_R = 0, dm_R = 8, m_R = 8, initfreq_R = c(1000, 1700),
               gmode = gmode,
               thruth_data=true_data, actPlot=TRUE,
               limFreq = c(1000));
    
    result[i+(j-1)*N,1]=out$covpbb[1,1]
    result[i+(j-1)*N,2]=dist

    result[i+(j-1)*N,3]=out$covpbb[1,2]
    result[i+(j-1)*N,4]=out$covpbb[1,3]
    result[i+(j-1)*N,5]=out$residual[1,2]
    result[i+(j-1)*N,6]=out$residual[1,3]
    result[i+(j-1)*N,7]=out$residual[1,4]

    result[i+(j-1)*N,8]=out1$covpbb[1,2]
    result[i+(j-1)*N,9]=out1$covpbb[1,3]
    result[i+(j-1)*N,10]=out1$residual[1,2]
    result[i+(j-1)*N,11]=out1$residual[1,3]
    result[i+(j-1)*N,12]=out1$residual[1,4]
    print(out1$covpbb)
  }
}

print(out)
print(out1)

filename=sprintf("results_%s_fcut1000Hz_AA_%s.txt", filtering_method, name)
write.table(result, file=filename, sep=" ", row.names=FALSE, col.names=FALSE)

#########################
### BAYESIAN ANALYSIS ###
#########################
# determine the time interval where signal is non null
#s20_0 = s20[,2]!=0;  # identifying values == 0
#s20_0 = s20[s20_0,]; # discarding values

#timeGmode = s20_0[c(1,length(s20_0[,1])), 1];#time interval for signal different from 0  
# postSamples = read.table("postsamples.txt", header = T); # importing posterios samples
# 
# covpbbBayes(noisydata, postSamples, l=200, p=90, fs=fs, movGmode = 11, 
#             um = 10, dm = 10, movBand = 5, timeGmode = NULL, 
#             Obergaunlinguer_data, actPlot =TRUE, fmean = mean(fits_data$f))


