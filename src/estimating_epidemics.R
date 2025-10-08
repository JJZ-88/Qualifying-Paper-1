library(rtestim)
library(ggplot2)
library(signal)


##Simulated epidemic scenarios
get_R_true <- function(scen_num, n_days){
  R_true <- c() #true instantaneous reproduction numbers
  if(scen_num == 1){ #piecewise constant
    R_true <- c(rep(2.5, 70), rep(0.5, 160), rep(2.5, n_days - 230))
  }
  else if(scen_num == 2){ #sinusoidal
    R_true <- 1.3 + 1.2 * sin(pi * seq(1, n_days, by = 1) / 60)
  }
  else if(scen_num == 3){ #exponential rise and fall
    R_true <- c(exp(0.03 * seq(1, 40, by = 1)), exp(-0.015 * seq(1, 150, by = 1) + 1.2), 
                exp(0.02 * seq(1, n_days - 190, by = 1) - 1.05))
  }
  
  return(R_true)#list('R_true' = R_true, 'scen' = scen_num))
}

#get 'probability mass' of gamma distribution (serial time dist.)
get_dist_vals <- function(dist_param, n_days){
  dist <- pgamma(seq(1, n_days, by = 1), dist_param[1], rate = dist_param[2]) - 
    pgamma(seq(0, n_days - 1, by = 1), dist_param[1], rate = dist_param[2])
  return(rev(dist)) #return dist in reverse time, so that w_1 is in the last index
}

#get the true R_t, incidence, infectiousness numbers
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
      #if(i < 100){
      #  print(dist[(n_days - i + 2):n_days])
      #}
      
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


n_days <- 300
#dist_param <- c((1 / 0.65)^2, (1 / 0.65)^2 / 6.5)
dist_param <- c(2.7066, 2.7066 / 15.3)
init_inc <- 10
offset <- 20

epi_vals <- setup(2, n_days, dist_param, init_inc, offset) #get true values
epi_vals
r_lam <- diff(log(epi_vals$lam_true)) #r_t|lam_t as pairwise difference of log(lam_t)
tau <- 8
#get the true r_t from true R_t based on same model assumptions
r_true <- (epi_vals$R_true^(1 / dist_param[1]) - 1) * dist_param[2]

rtestim_lambda <- cv_estimate_rt(epi_vals$I_true + 1, nfold = 5, korder = 3)
rtestim_fit <- fitted(rtestim_lambda, 'lambda.min')
#rtestim_fit <- estimate_rt(epi_vals$I_true + 0.001, lambda = 0)

plot(seq(1, 280, 1), epi_vals$R_true, type = 'l')
lines(rtestim_fit)

I_emp <- sgolayfilt(epi_vals$I_true, p = 3, n = 51)
r_emp <- sgolayfilt(diff(log(I_emp)), p = 3, n = 91)

#plot 2: growth rates
plot(seq(1, 280, 1), r_true, type = 'l')
lines(seq(tau + 1, 280, 1), r_emp[1:(length(r_emp) - tau + 1)], col = 'blue')
lines(seq(1, 280 - tau, 1), r_lam[tau:length(r_lam)], col = 'red')
lines(seq(1, 280, 1), 
      (rtestim_fit^(1 / dist_param[1]) - 1) * dist_param[2])




  
  
  
  
  
