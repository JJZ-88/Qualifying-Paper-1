library(rtestim)
library(ggplot2)
library(signal)
library("gridExtra")
library(ggpubr)
library(latex2exp)


##Simulated epidemic scenarios
get_R_true <- function(scen_num, n_days){
  R_true <- c() #true instantaneous reproduction numbers
  if(scen_num == 1){ #piecewise constant
    R_true <- c(rep(2, 100), rep(0.7, 100), rep(1.3, n_days - 200))
  }
  else if(scen_num == 2){ #sinusoidal
    R_true <- 1.3 + 1.2 * sin(pi * seq(1, n_days, by = 1) / 60)
  }
  else if(scen_num == 3){ #exponential rise and fall
    R_true <- c(exp(0.03 * seq(1, 40, by = 1)), exp(-0.015 * seq(1, 150, by = 1) + 1.2), 
                exp(0.02 * seq(1, n_days - 190, by = 1) - 1.05))
  }
  else if(scen_num == 4){
    #R_true <- c(1.5 + 0.01 * seq(1, 150, 1), 3 - 0.01 * seq(1, 150, 1))
    R_true <- c(seq(1, 100, 1) / 200 + sin(pi * seq(1, 100, 1) / 20 )/10 + 1.5, 
                rep(0.8, 100), exp(0.005 * seq(1, 100, 1)) * 0.8)
  }
  return(R_true)#list('R_true' = R_true, 'scen' = scen_num))
}


##'Probability Mass' of gamma distribution (serial time dist.)
get_dist_vals <- function(dist_param, n_days){
  dist <- pgamma(seq(1, n_days, by = 1), dist_param[1], rate = dist_param[2]) - 
    pgamma(seq(0, n_days - 1, by = 1), dist_param[1], rate = dist_param[2])
  return(rev(dist)) #return dist in reverse time, so that w_1 is in the last index
}

##Get true R_t, incidence, infectiousness numbers
setup <- function(scen_num, n_days, dist_param, init_inc, offset){
  R_true <- get_R_true(scen_num, n_days) #get R_t and w_t based on preset parameters
  dist <- get_dist_vals(dist_param, n_days)
  I_true <- lam_true <- rep(0, n_days)
  I_true[1] <- init_inc #initialize initial incidence
  small_epi <- T
  
  #loop until the epidemic is sufficiently 'large', i.e. total incidences >= 500
  while(small_epi == T){ 
    for(i in 2:n_days){ #iterate over days
      #Lam_t = sum_{j=1}^{i-1} I_{i-j} * w_j 
      lam_true[i] <- sum(I_true[1:i - 1] * dist[(n_days - i + 2):n_days])
      
      #I_t ~ Poisson(lam_t * R_t)
      I_true[i] <- rpois(1, lam_true[i] * R_true[i])
    }
    if(sum(I_true) >= 500){
      small_epi = F
    }
  }
  return(list('R_true' = R_true[(offset+1):n_days], 'I_true' = I_true[(offset+1):n_days], 
              'lam_true' = lam_true[(offset+1):n_days]))
}

##Run rtestim package to get optimal (tuned) fit
run_rtestim <- function(incidence, dist_param, order){
  rtestim_lambda <- cv_estimate_rt(incidence, nfold = 5, korder = order, dist_gamma = c(dist_param[1], 1 / dist_param[2]))
  rtestim_fit <- fitted(rtestim_lambda, 'lambda.min')
  rtestim_ci <- confband(rtestim_lambda, lambda = 'lambda.min')
  return(list('lambda' = rtestim_lambda$lambda.min, 'fit' = rtestim_fit, 'ci' = rtestim_ci))
}
test_for_daniel <- simulate_epidemic(n_days, dist_param, init_inc, 0, offset, 1, 0, "")

##Compute r_t from R_t from Wallinga Lipsitch
get_rt <- function(R, dist_param){
  return((R^(1 / dist_param[1]) - 1) * dist_param[2])
}

simulate_epidemic <- function(n_days, dist_param, init_inc, order, offset, scen_num, save_plot, file_name){
  ##Derive true values and model estimates
  epi_vals <- setup(scen_num, n_days, dist_param, init_inc, offset) #get true values
  r_lam <- diff(log(epi_vals$lam_true)) #r_t | lam_t as pairwise difference of log(lam_t)
  tau <- round(dist_param[1] / dist_param[2] / 2)
  #get the true r_t from true R_t based on same model assumptions
  r_true <- get_rt(epi_vals$R_true, dist_param)
  #esimate Rt using rtestim package
  R_est <- run_rtestim(epi_vals$I_true, dist_param, order)
  #smoothed incidence using SG filter (window 50, cubic fit)
  I_emp <- pmax(sgolayfilt(epi_vals$I_true, p = 3, n = 31), rep(0.00001, n_days - offset)) 
  r_emp <- sgolayfilt(diff(log(I_emp)), p = 3, n = 91) #r_t | S_t from smoothed incidence
  r_Rt <- get_rt(R_est$fit, dist_param)
  r_upperCI <- get_rt(R_est$ci$`97.5%`, dist_param)
  r_lowerCI <- get_rt(R_est$ci$`2.5%`, dist_param)
  
  
  ##Plot 1: estimating reproduction numbers
  g1 <- ggplot(data.frame(x = seq(offset + 1, n_days, 1), y = epi_vals$R_true, z = R_est$fit)) +
    geom_line(aes(x = x, y = y, colour = 'R_true')) +
    geom_line(aes(x = x, y = z, colour = 'R_estimate')) +
    #geom_line(aes(x = x, y = R_est$ci$`2.5%`, colour = 'R_estimate'), linetype = 'dotted') +
    #geom_line(aes(x = x, y = R_est$ci$`97.5%`, colour = 'R_estimate'), linetype = 'dotted') +
    scale_color_manual(values = c('R_true' = 'black', 'R_estimate' = 'blue')) +
    theme(legend.position = "bottom") +
    labs(title = 'Estimating Instantaneous \n Reproduction Number', x = 'Days', y = 'Reproduction Number',
         color = '')

  
  ##Plot 2: estimating growth rates 3 ways
  g2 <- ggplot() + 
    geom_line(aes(x = seq(offset + 1, n_days, 1), y = r_true, colour = 'r_true')) +
    geom_line(aes(x = seq(offset + tau + 1, n_days, 1), y = r_emp[1:(length(r_emp) - tau + 1)], colour = 'r_t | S_t')) +
    geom_line(aes(x = seq(offset + 1, n_days - tau, 1), y = r_lam[tau:length(r_lam)], colour = 'r_t | lam_t')) +
    geom_line(aes(x = seq(offset + 1, n_days, 1), y = r_Rt, colour = 'r_t | R_t')) +
    #geom_line(aes(x = seq(offset + 1, n_days, 1), y = r_upperCI, colour = 'r_t | R_t'), linetype = 'dotted') +
    #geom_line(aes(x = seq(offset + 1, n_days, 1), y = r_lowerCI, colour = 'r_t | R_t'), linetype = 'dotted') +
    scale_color_manual(values = c('r_true' = 'black', 'r_t | S_t' = 'blue', 'r_t | lam_t' = 'red', 'r_t | R_t' = 'green')) +
    theme(legend.position = "bottom") +
    labs(title = 'Estimating Instantaneous \n Growth Rate',x = 'Days', y = 'Growth Rate', 
         color = '')
  
  ##Plot 3: incidence
  #Here, we use I_t coming from a Poisson distribution with rate (R_t lam_t)
  g3 <- ggplot(data.frame(x = seq(offset + 1, n_days, 1), y = epi_vals$I_true, z = epi_vals$lam_true * R_est$fit)) + 
    geom_point(aes(x = x, y = y, colour = 'I_true'), size = 0.2) +
    geom_line(aes(x = x, y = z, colour = 'I_t | R_t')) +
    scale_color_manual(values = c('I_true' = 'black', 'I_t | R_t' = 'blue')) + 
    theme(legend.position = "bottom") +
    labs(x = 'Days', y = 'Incidence of Infection', title = 'Estimating Incidence \n of Infection', 
         color = '')
  
  ##Plot 4: smoothed incidence curves
  g4 <- ggplot(data.frame(x = seq(offset + 1, n_days - 2 * tau, 1), y = epi_vals$I_true[1:(n_days - offset - 2 * tau)], 
                    z = I_emp[1:(n_days - offset - 2 * tau)], w = epi_vals$lam_true[(2 * tau + 1):(n_days - offset)])) + 
    geom_point(aes(x = x, y = y, colour = 'I_t'), size = 0.2) +
    geom_line(aes(x = x, y = z, colour = 'S_t')) +
    geom_line(aes(x = x, y = w, colour = 'lam_{t + 2tau}')) +
    scale_color_manual(values = c('I_t' = 'black', 'S_t' = 'blue', 'lam_{t + 2tau}' = 'red')) + 
    theme(legend.position = "bottom") +
    labs(x = 'Days', y = 'Incidence of Infection', title = 'Smoothing Incidence \n of Infection', 
         color = '')
  
  print(g1)
  print(g2)
  print(g3)
  print(g4)
  
  if(save_plot == 1){
    to_plot <- ggarrange(g1, g2, g3, g4, nrow = 2, labels = c('a)', 'b)', 'c)', 'd)'))
    ggsave(file_name, to_plot)
  }
  if(save_plot == 2){
    to_plot <- ggarrange(g1, g2, nrow = 1, labels = c('a)', 'b)'))
    ggsave(file_name, to_plot)
  }
  
  return(list('R_true' = epi_vals$R_true, 'I_true' = epi_vals$I_true, 'r_est' = r_Rt,
              'lam_true' = epi_vals$lam_true,'r_true' = r_true, 'R_est' = R_est$fit, 
              'g1' = g1, 'g2' = g2, 'g3' = g3, 'g4' = g4))
}

kl <- function(R_true, R_est, lambda_true){
  return(sum(lambda_true * (R_true * log(R_true / R_est) + R_est - R_true)))
}

mse <- function(r_true, r_est){
  return(sum((r_true - r_est)^2))
}


##Initialize parameters
set.seed(88)
n_days <- 300
#dist_param <- c((1 / 0.65)^2, (1 / 0.65)^2 / 6.5)
dist_param <- c(2.7066, 2.7066 / 15.3)
init_inc <- 10
offset <- 20
order <- 3
scen_num <- 2

setwd("/Users/JZ/Desktop/Qualifying-Paper-1/report/fig")
epi1 <- simulate_epidemic(n_days, dist_param, init_inc, order, offset, scen_num, 0, 'epi_paper.pdf')

kl(epi1$R_true, epi1$R_est, epi1$lam_true)
mse(epi1$r_true, epi1$r_est)


##Estimates with misspecified serial interval distribution
dist_param_misspec <- c(2.7066, 2.7066 / 15.3 * 3) #misspecify gamma mean as 1/3 of true mean
dist_misspec <- get_dist_vals(dist_param_misspec, n_days)

Rt_misspec <- run_rtestim(epi1$I_true, dist_param_misspec, order)

kl(epi1$R_true, Rt_misspec$fit, epi1$lam_true)
mse(epi1$r_true, get_rt(Rt_misspec$fit, dist_param_misspec))

##Plot 5: estimating reproduction numbers under misspecification of generation time distribution
irr_misspec <- ggplot(data.frame(x = seq(offset + 1, n_days, 1), y = epi1$R_true, z = epi1$R_est, w = Rt_misspec$fit)) +
  geom_line(aes(x = x, y = y, colour = 'R_true')) +
  geom_line(aes(x = x, y = z, colour = 'R_est')) +
  geom_line(aes(x = x, y = w, colour = 'R_est (misspecified)')) +
  geom_line(aes(x = x, y = Rt_misspec$ci$`2.5%`, colour = 'R_est (misspecified)'), linetype = 'dotted') +
  geom_line(aes(x = x, y = Rt_misspec$ci$`97.5%`, colour = 'R_est (misspecified)'), linetype = 'dotted') +
  scale_color_manual(values = c('R_true' = 'black', 'R_est' = 'blue', 'R_est (misspecified)' = 'red')) + 
  theme(legend.position = "bottom") +
  labs(title = 'Estimating Instantaneous Reproduction Number \n under Misspecied GTD', 
       x = 'Days', y = 'Reproduction Number', color = '', tag = 'a')

##Plot 6: estimating growth rates under misspecification of serial interval distribution
gr_misspec <- ggplot(data.frame(x = seq(offset + 1, n_days, 1), y = epi1$r_true, z = get_rt(Rt_misspec$fit, dist_param_misspec), 
                  w = get_rt(Rt_misspec$ci$`97.5%`, dist_param_misspec), v = get_rt(Rt_misspec$ci$`2.5%`, dist_param_misspec))) + 
  geom_line(aes(x = x, y = y, colour = 'r_true')) +
  geom_line(aes(x = x, y = z, colour = 'r_t | R_t')) +
  geom_line(aes(x = x, y = w, colour = 'r_t | R_t'), linetype = 'dotted') +
  geom_line(aes(x = x, y = v, colour = 'r_t | R_t'), linetype = 'dotted') +
  scale_color_manual(values = c('r_true' = 'black', 'r_t | R_t' = 'blue')) +
  theme(legend.position = "bottom") +
  labs(title = 'Estimating Instantaneous Growth Rate',x = 'Days', y = 'Growth Rate', color = '', tag = 'b')

to_plot_misspec <- ggarrange(irr_misspec, gr_misspec, ncol = 2, labels = c('a)', 'b)'))
ggsave('epi_misspec.pdf', to_plot_misspec)

##Plot Piecewise Constant Epidemic
pc_0 <- simulate_epidemic(n_days, dist_param, init_inc, 0, offset, 1, 0, "pc_0.pdf")
pc_3 <- simulate_epidemic(n_days, dist_param, init_inc, 3, offset, 1, 0, "pc_3.pdf")

to_plot_pc <- ggarrange(pc_0$g1, pc_0$g2, pc_3$g1, pc_3$g2, ncol = 2, nrow = 2, labels = c('a)', 'b)', 'c)', 'd)'))
ggsave('epi_pc.pdf', to_plot_pc)

kl(pc_0$R_true, pc_0$R_est, pc_0$lam_true)
kl(pc_3$R_true, pc_3$R_est, pc_3$lam_true)
mse(pc_0$r_true, pc_0$r_est)
mse(pc_3$r_true, pc_3$r_est)

##Plot Composite Epidemic
comp_1 <- simulate_epidemic(n_days, dist_param, init_inc, 1, offset, 4, 0, "comp_0.pdf")
comp_3 <- simulate_epidemic(n_days, dist_param, init_inc, 3, offset, 4, 0, "comp_3.pdf")

to_plot_comp <- ggarrange(comp_1$g1, comp_1$g2, comp_3$g1, comp_3$g2, ncol = 2, nrow = 2, labels = c('a)', 'b)', 'c)', 'd)'))
ggsave('epi_comp.pdf', to_plot_comp)

kl(comp_1$R_true, comp_1$R_est, comp_1$lam_true)
kl(comp_3$R_true, comp_3$R_est, comp_3$lam_true)


##Penalized r_t - mini proposal

#make nth differences matrix
make_div_diff <- function(n, it){
  d <- cbind(rep(0, n - 1), diag(rep(1, n - 1))) + cbind(diag(rep(-1, n - 1)), rep(0, n - 1)) #first differences
  if(it > 1){ 
    for(i in 1:(it - 1)){ #for nth differences greater than 2, make 1st differences matrix of 1 dimension lower and multiply to d
      d0 <- cbind(rep(0, n - 1 - i), diag(rep(1, n - 1 - i))) + cbind(diag(rep(-1, n - 1 - i)), rep(0, n - 1 - i))
      d <- d0 %*% d
    }
  }
  return(d)
}

#objective function to minimize
objective_fn <- function(r, inc, lambda_1, lambda_2, d, d2){
  n <- length(inc)
  loss <- sum(-inc[2:n]*(log(inc[1:(n - 1)]) + r[2:n]) + inc[1:(n - 1)] * exp(r[2:n])) #poisson loss
  pen <- lambda_1 * sum(exp((d %*% r)^2))  #penalize first differences
  pen_2 <- lambda_2 * sum(exp((d2 %*% r)^2)) #penalize second differences
  #pen <- lambda_1 * sum((exp(r[2:n]) - exp(r[1:(n-1)]))^2)
  #pen <- lambda * sum((d %*% r)^2)
  return(loss + pen + pen_2)
}

#use optim function to minimize objective
#NOTE EACH CALL TO OBJECTIVE NEEDS TO CHANGE FUNCTION objective_fn ACCORDINGLY
obj_fs = optim(rep(0, 280), fn = objective_fn, inc = e, d = make_div_diff(280, 1), d2 = make_div_diff(280, 2), lambda_1 = 200, lambda_2 = 500, method = 'BFGS')

plotr_fsdiff <- ggplot(data.frame(x = seq(tau + 1, n_days - offset, 1), y = obj_fs$par[1:(n_days - offset - tau)],z = epi1$r_true[(tau + 1):(n_days - offset)])) +
  geom_line(aes(x = x, y = y, colour = 'r_est')) +
  geom_line(aes(x = x, y = z, colour = 'r_true')) +
  scale_color_manual(values = c('r_true' = 'black', 'r_est' = 'blue')) + 
  theme(legend.position = "bottom") +
  labs(title = 'Smoothed Estimator of Growth Rate \n
       with First and Second Differences',x = 'Days', y = 'Growth Rate', color = '', tag = 'c')

plot(plotr_fsdiff)
obj_f = optim(rep(0, 280), fn = objective_fn, inc = e, d = make_div_diff(280, 1), d2 = make_div_diff(280, 2), lambda_1 = 200, lambda_2 = 500, method = 'BFGS')

plotr_fdiff <- ggplot(data.frame(x = seq(tau + 1, n_days - offset, 1), y = obj_f$par[1:(n_days - offset - tau)],z = epi1$r_true[(tau + 1):(n_days - offset)])) +
  geom_line(aes(x = x, y = y, colour = 'r_est')) +
  geom_line(aes(x = x, y = z, colour = 'r_true')) +
  scale_color_manual(values = c('r_true' = 'black', 'r_est' = 'blue')) + 
  theme(legend.position = "bottom") +
  labs(title = 'Smoothed Estimator of Growth Rate \n
       with First Differences',x = 'Days', y = 'Growth Rate', color = '', tag = 'b')

obj_mle = optim(rep(0, 280), fn = objective_fn, inc = e, d = make_div_diff(280, 1), d2 = make_div_diff(280, 2), lambda_1 = 200, lambda_2 = 500, method = 'BFGS')

plotr_mle <- ggplot(data.frame(x = seq(tau + 1, n_days - offset, 1), y = obj_mle$par[1:(n_days - offset - tau)],z = epi1$r_true[(tau + 1):(n_days - offset)])) +
  geom_line(aes(x = x, y = y, colour = 'r_est')) +
  geom_line(aes(x = x, y = z, colour = 'r_true')) +
  scale_color_manual(values = c('r_true' = 'black', 'r_est' = 'blue')) + 
  theme(legend.position = "bottom") +
  labs(title = 'Smoothed Estimator of Growth Rate \n
       with MLE',x = 'Days', y = 'Growth Rate', color = '', tag = 'a')
  
to_plot_smoothrt <- grid.arrange(plotr_mle, plotr_fdiff, plotr_fsdiff, nrow = 2)
ggsave('estimate_growthrate.pdf', to_plot_smoothrt)

#used for basic tuning of lambda parameters
lambda_test <- seq(400, 800, by = 20)
lambda_mse <- seq(400, 800, by = 20)
for(l in 1:length(lambda_test)){
  obj = optim(rep(0, 280), fn = f, inc = e, d = make_div_diff(280, 1), d2 = make_div_diff(280, 2), lambda = lambda_test[l], method = 'BFGS')
  lambda_mse[l] <- mse(obj$par[1:(n_days - offset - tau)], epi1$r_true[(tau + 1):(n_days - offset)])
}
data.frame(x = lambda_test, mse = lambda_mse)





