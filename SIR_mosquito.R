library(ggplot2)
library(deSolve)
library(tidyverse)
library("outbreaks")

# SIR model with population structure and vector transmission
sir_model_mosquito <- function(times, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    N = S1 + S2 + I1 + I2 + R1 + R2 # Total population size
    
    # Transmission rates depend on the number of mosquitoes
    beta11_mosq <- beta11 * N_m
    beta12_mosq <- beta12 * N_m
    beta21_mosq <- beta21 * N_m
    beta22_mosq <- beta22 * N_m
    
    dS1 <- -beta11_mosq * S1 * I1 - beta12_mosq * S1 * I2
    dS2 <- -beta22_mosq * S2 * I2 - beta21_mosq * S2 * I1
    dI1 <-  beta11_mosq * S1 * I1 + beta12_mosq * S1 * I2 - gamma * I1
    dI2 <-  beta22_mosq * S2 * I2 + beta21_mosq * S2 * I1 - gamma * I2
    dR1 <-                                        gamma * I1
    dR2 <-                                        gamma * I2
    
    return(list(c(dS1, dS2, dI1, dI2, dR1, dR2)))
  })
}


# set initial conditions, parameters,
# sequence of time points for solving the model

init = c(S1 = 3695, S2 = 3694, I1 = 1, I2 = 1, R1 = 0, R2 = 0) 
# split population into half for initial susceptibles?

parameters = c(beta11 = 0.5/9725,
               beta12 = 0.025/9725,
               beta21 = 0.025/9725,
               beta22 = 1/9725,
               gamma = 1/5,
               N_m = 200)

times = seq(0, 200, length.out=1000)

out = ode(y = init, times = times, func = sir_model, parms = parameters)
out = as.data.frame(out)

# plot both SIR models on the same graph
out_long = pivot_longer(out, cols=c(2,3,4,5,6,7))
ggplot(out_long, aes(x=time, y=value, group=name, colour=name)) +
  geom_line() +
  theme_minimal() +
  xlab("Days Elapsed") +
  ylab("Number of individuals") +
  guides(colour=guide_legend(title="Compartment")) +
  ggtitle("SIR model mosquito with homogeneity")

# We are interested in the TOTAL number of infecteds
# For each time, we should add the values for I1 and I2 in out
out$total_I = out$I1 + out$I2

out_long = pivot_longer(out, cols=c(4,5))
ggplot() +
  geom_line(data=out_long, aes(x=time, y=total_I, group=name, colour='Simulated')) +
  geom_line(data=zika_yap_2007, aes(x=nr, y=value, color='Real Data')) +
  theme_minimal() +
  xlab("Days Elapsed") +
  ylab("Number of individuals") +
  ggtitle("SIR Model With Population Structure: Cases and vector transmission")
