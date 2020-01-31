# Checking different linear models to explain general relationshiop
# between frequency and ratios

#############
### Model ###
#############

fits_data = read.table("fits_data.dat", sep = ",");# data to generate model
colnames(fits_data) = c("f", "r");

#mod = lm(r~ 0 + f + I(f^2) + I(f^3), data = fits_data); # first model
#summary(mod);
mod = lm(r~ 0 + f + I(f^3), data = fits_data);# final model
summary(mod);
AIC(mod);BIC(mod);

res = fit$residuals
fit = lm (r~ 0 + f + I(f^2) + I(f^3), data=fits_data, weights=(1/fits_data$f^2))
summary(fit)
plot(fit)

# This model allows r>0
fit = lm(log(r)~ 0 + f + I(f^2) + I(f^3), data = fits_data); # first model
summary(fit)
plot(fit)

fit = lm(r~ 0 + log(f) + I(log(f)^2) + I(log(f)^3), data = fits_data); # first model
summary(fit)
plot(fit)

fit = lm(log(r)~ 0 + log(f) + I(log(f)^2) + I(log(f)^3), data = fits_data); # first model
summary(fit)
fit = lm(log(r)~ 0 + log(f) + I(log(f)^2), data = fits_data); # first model
summary(fit)
plot(fit)


##############################
### Box-Cox Transformation ###
##############################

# https://datascienceplus.com/how-to-detect-heteroscedasticity-and-rectify-it/

library(caret)
bc_mod    = caret::BoxCoxTrans(fits_data$r);
fits_data = cbind(fits_data, bc_r = predict(bc_mod, fits_data$r));

bc_mod = lm(bc_r ~ f, data=fits_data);
plot(bc_mod);
summary(bc_mod);

lmMod <- lm(bc_r ~ 0 + f + I(f^3), data=fits_data)
plot(bc_mod);
summary(bc_mod);

#############
### lmvar ###
#############

library(lmvar)

mu1 = "~ f - 1";
mu2 = "~ f + I(f^2) - 1";
mu3 = "~ f + I(f^2) + I(f^3) - 1";
mu4 = "~ f + I(f^3) - 1";

s1 = "~ f - 1"; 
s2 = "~ f + I(f^2) - 1";

mu = c(mu1, mu2, mu3, mu4);
s  = c(s1, s2);

aic = matrix(NA, nrow = length(s), ncol = length(mu));
bic = matrix(NA, nrow = length(s), ncol = length(mu));

count_i = count_j = 0;

for(i in s){
  
  count_i = count_i + 1;
  count_j = 0;
  
  for(j in mu){
    
    count_j = count_j + 1;
    
    Xm  = model.matrix(eval(parse(text=j)), fits_data);
    Xs  = model.matrix(eval(parse(text=i)), fits_data);
    
    fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = FALSE);
    
    aic[count_i, count_j] = AIC(fit);
    bic[count_i, count_j] = BIC(fit);
    
  }
}

aic
bic

### plot ###

Xm  = model.matrix(eval(parse(text=mu3)), fits_data);
Xs  = model.matrix(eval(parse(text=s2)), fits_data);
fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = FALSE);

plot(fits_data$f, fits_data$r);
pred = predict(fit, X_mu = Xm, X_sigma = Xs, interval = "prediction"); # predictions
lwr = cbind(fits_data$f, pred[,3]); # lower band
lwr = lwr[order(lwr[,1]), ];
lines(lwr[,1], lwr[,2], col = "red", lwd = 2);
upr = cbind(fits_data$f, pred[,4]); # upper band
upr = upr[order(upr[,1]), ];
lines(upr[,1], upr[,2], col = "red", lwd = 2);
med = cbind(fits_data$f, pred[,1]); # middle point
med = med[order(med[,1]), ];
lines(med[,1], med[,2], col = "blue", lwd = 2);

# Model with intercept
Xm  = model.matrix(eval(parse(text=mu3)), fits_data);
Xs  = model.matrix(eval(parse(text=s2)), fits_data);

fit = lmvar(fits_data$r, X_mu = Xm, X_sigma = Xs, intercept_mu = TRUE);

aic[count_i, count_j] = AIC(fit);
bic[count_i, count_j] = BIC(fit);