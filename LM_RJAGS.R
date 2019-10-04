library('rjags')

fits <- read.table(file="inputs/fits_data.dat", header=FALSE, sep=",")
f <- fits[,1] # maybe better to center f already here and then use uncentered jags code
r <- fits[,2]
n<-nrow(fits)


fnew<- seq(45,1700,length=10)
K<-length(fnew)
# to be substituted by the values of the g-modes

# Put the BUGS code in a character string
BUGSmodel = "model
{
  for (i in 1:n){
  mu[i] <- a0 + a1*(f[i]-mean(f))+a2*(f[i]-mean(f))^2+a3*(f[i]-mean(f))^3
 # mu[i] <- a1*f[i]+a2*f[i]^2+a3*f[i]^3
  sigma2[i] <-  b1*f[i] 
  tau[i] <- 1/sigma2[i]
  r[i] ~ dnorm(mu[i],tau[i]) 
  }
#prediction
  for (i in 1:K){
   munew[i] <- a0+ a1*(fnew[i]-mean(f))+a2*(fnew[i]-mean(f))^2+a3*(fnew[i]-mean(f))^3
  # munew[i] <- a1*fnew[i]+a2*fnew[i]^2+a3*fnew[i]^3
  sigma2new[i] <-  b1*fnew[i] 
  tau2new[i] <- 1/sigma2new[i]
  rnew[i] ~ dnorm(munew[i],tau2new[i])
  }
  a0 ~ dnorm(0,0.001)
  a1 ~ dnorm(0,0.001)
  a2 ~ dnorm(0,0.001)
  a3 ~ dnorm(0,0.001)
  #b0 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
}
"



data = list(n=n,r=r,f=f,fnew=fnew,K=K)
# The initial values as a list (JAGS will automatically generate these if unspecified)
#inits = list(a1=1.0e-6, a2=1.0e-9, a3=1.0e-12, b1=1.0e-10)
inits = list(a0=0,a1=-1.0e-6, a2=6.0e10, a3=1.0e-13, b1=1.0e-10) # for mean-centered version
# parameters to monitor
parameters = c('a0','a1','a2', 'a3', 'b1','rnew')

# How many burn-in steps?
burn_in = 1000

# How many proper steps?
steps = 10000

# Thinning?
thin = 1

# Random number seed
seed = 42



#compilation of the BUGS model, no Gibbs sampling yet
foo <- jags.model(textConnection(BUGSmodel),data=data, inits=inits) 

#burnin samples
update(foo, burn_in)

#draws n.iter MCMC samples, monitors the parameters specified in variable.names, thins the
#output by a factor of thin and stores everything in a MCMC.list object
out <- coda.samples(model=foo, variable.names=parameters, n.iter=steps, thin=thin)





#Obtain summary statistics for each of the monitored parameters
summary(out)

#Obtain traceplots and kernel density estimates for each parameter
par(mar = rep(2, 4)) #adjusts margin plotting parameters
plot(out)


outmatrix <- as.matrix(out)

a0hat <- median(outmatrix[,1])
a1hat <- median(outmatrix[,2])
a2hat <- median(outmatrix[,3])
a3hat <- median(outmatrix[,4])
rhat <- a0hat + a1hat*(f-mean(f))+a2hat*(f-mean(f))^2+ a3hat*(f-mean(f))^3
#rhat <- a1hat*f+a2hat*f^2+ a3hat*f^3
#rnewhat <- a1hat*fnew + a2hat*fnew^2 + a3hat*fnew^3
rnewhat <- a0hat + a1hat*(fnew-mean(f)) + a2hat*(fnew-mean(f))^2 + a3hat*(fnew-mean(f))^3

par(mfrow = c(1,1))
plot(f, rhat)
points(f, r, col='red')

plot(fnew, rnewhat, col='black', type='l')
points(f, r, col='red')

write.table(outmatrix, file="postsamples.txt", row.names=FALSE, col.names=TRUE)
