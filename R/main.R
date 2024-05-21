library(arrow) 
library(rcpm)
library(marqLevAlg)
library(ggplot2)
library(tibble)

devtools::load_all('C:/Users/Yusuf/Documents/rcpm')

#As the rcpm package by Corentin Segalas is not yet on CRAN we hope to troubleshoot locally by reloading after modifying the code. 
castem <- read_feather('data/nested_case_control.feather')
esti_fixedclass <- rcpme(castem, isa15 ~ delaidiag | ID2, nbnodes = 10, model = "isplines", statut = "statut")
esti_latentclass1 <- rcpme(castem, isa15 ~ delaidiag | ID2, nbnodes = 10, model = "isplines",link = "linear", statut = "statut", latent = TRUE, classprob = ~ 1, lambda = 100)
save(esti_latentclass1, file = sprintf("C:/Users/Yusuf/Documents/rcpm/models/opt_lambda_%g.RData", lambda))
save(esti_fixedclass, file = "trained models/esti_fixedclass.RData") #
save(esti_latentclass1, file = "trained models/esti_latentclass1.RData") #

esti_latentclassCEP <- rcpme(castem, isa15 ~ delaidiag | ID2, nbnodes = 10, model = "isplines", statut = "statut", classprob = ~ CEP)


# 'esti_fixedclass'
vis <- seq(-25,0,by=0.1)
ispl = matrix(rep(NA,3*length(vis)), nrow = length(vis), ncol = 3)
for (i in seq(length(vis))) {
  ispl[i,] <- ispline(vis[i] - esti_fixedclass$optpar[3], 0, 20, 5)}
dat <- cbind(vis, esti_fixedclass$optpar[1] + esti_fixedclass$optpar[2]*vis, esti_fixedclass$optpar[1] + esti_fixedclass$optpar[2]*vis - ispl %*% esti_fixedclass$optpar[12:14]^2)
colnames(dat) = c("x","y","z")
dat <- as_tibble(dat)
ggplot(data = dat) + geom_line(aes(x = x, y = y), colour = "blue") + geom_line(aes(x=x,y=z), colour = "red")


# 'esti_latentclass1'
# prob
exp(esti_latentclass1$optpar[15])/(1+exp(esti_latentclass1$optpar[15]))
# courbes marginales
vis <- seq(-25,0,by=0.1)
ispl = matrix(rep(NA,3*length(vis)), nrow = length(vis), ncol = 3)
for (i in seq(length(vis))) {
  ispl[i,] <- ispline(vis[i] - esti_latentclass1$optpar[3], 0, 20, 5)}
dat <- cbind(vis, esti_latentclass1$optpar[1] - esti_latentclass1$optpar[2]*vis, esti_latentclass1$optpar[1] - esti_latentclass1$optpar[2]*vis - ispl %*% esti_latentclass1$optpar[12:14]^2)
colnames(dat) = c("x","y","z")
dat <- as_tibble(dat)
ggplot(data = dat) + geom_line(aes(x = x, y = y), colour = "blue") + geom_line(aes(x=x,y=z), colour = "red")

# 'esti_latentclassCEP'
# prob
exp(esti_latentclassCEP$optpar[15])/(1+exp(esti_latentclassCEP$optpar[15]))
# courbes marginales
vis <- seq(-25,0,by=0.1)
ispl = matrix(rep(NA,3*length(vis)), nrow = length(vis), ncol = 3)
for (i in seq(length(vis))) {
  ispl[i,] <- ispline(vis[i] - esti_latentclassCEP$optpar[3], 0, 20, 5)}
dat <- cbind(vis, esti_latentclassCEP$optpar[1] - esti_latentclassCEP$optpar[2]*vis, esti_latentclassCEP$optpar[1] - esti_latentclassCEP$optpar[2]*vis - ispl %*% esti_latentclassCEP$optpar[12:14]^2)
colnames(dat) = c("x","y","z")
dat <- as_tibble(dat)
ggplot(data = dat) + geom_line(aes(x = x, y = y), colour = "blue") + geom_line(aes(x=x,y=z), colour = "red")


lambdas <- c(250, 300, 350, 400)
for (lambda in lambdas){
  print(lambda)
  fitted_model <- rcpme(castem, isa15 ~ delaidiag | ID2, nbnodes = 10, model = "linear-linear",link = "linear", statut = "statut", latent = TRUE, classprob = ~ 1, lambda = lambda)
  save(fitted_model, file = sprintf("C:/Users/Yusuf/Documents/latent-class-random-change-point-model/trained models/opt_lambda_%g.RData", lambda))
  print('appartition est:')
  print(exp(fitted_model$optpar[15]) / (1 + exp(fitted_model$optpar[15])))
  
}


spline <- read_feather('data/simulation_data_spline.feather')
lambdas_2 <- c(0,50,100) 

for(lambda in lambdas_2) {
  fitted_model <- rcpme(spline, y_trajectory ~ delaidiag | ID2, nbnodes = 10, model = "linear-linear",link = "linear", statut = "statut", latent = TRUE, classprob = ~ 1, lambda = lambda)
  save(fitted_model, file = sprintf("C:/Users/Yusuf/Documents/latent-class-random-change-point-model/trained models/new_simul_opt_lambda_%g.RData", lambda))
  print('appartition est:')
  print(exp(fitted_model$optpar[15]) / (1 + exp(fitted_model$optpar[15])))
}


load("C:/Users/Yusuf/Documents/latent-class-random-change-point-model/trained models/simul_opt_lambda_50.RData")
#viallonv@iarc.who.int
fitted_model
#Samuelson et al. 2007

library(dplyr)

summary_by_individual <- spline %>%
  group_by(ID2) %>%
  summarise(
    c_i = first(c_i),  # Assuming education remains constant
    beta_0i = mean(beta_0i),
    beta_2i = mean(beta_2i), 
    diff = first(age_at_diagnosis) - first(age)
  )

summary(summary_by_individual)