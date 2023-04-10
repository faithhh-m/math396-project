library(ggplot2)
library(deSolve)
library(tidyverse)
library("outbreaks")

# set up SEIR model

seir_model <- function(times, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    N = S1 + S2 + E1 + E2 + I1 + I2 + R1 + R2
    dS1 <- -beta11 * S1 * I1 - beta12 * S1 * I2
    dS2 <- -beta22 * S2 * I2 - beta21 * S2 * I1
    dE1 <-  beta11 * S1 * I1 + beta12 * S1 * I2 - sigma * E1
    dE2 <-  beta22 * S2 * I2 + beta21 * S2 * I1 - sigma * E2
    dI1 <-                                        sigma * E1 - gamma * I1
    dI2 <-                                        sigma * E2 - gamma * I2
    dR1 <-                                                     gamma * I1
    dR2 <-                                                     gamma * I2
    
    return(list(c(dS1, dS2, dE1, dE2, dI1, dI2, dR1, dR2)))
  })
}

# set initial conditions, parameters,
# sequence of time points for solving the model

init = c(S1 = 3695, S2 = 3694, E1 = 0, E2 = 0, I1 = 1, I2 = 1, R1 = 0, R2 = 0) 
# split population into half for initial susceptibles?

parameters = c(beta11 = 7/31120,
               beta12 = 7/38900,
               beta21 = 7/38900,
               beta22 = 7/15560,
               gamma = 1.43,
               sigma = 1/7)

times = seq(0, 200, length.out=1000)

out = ode(y = init, times = times, func = seir_model, parms = parameters)
out = as.data.frame(out)

# plot both SIR models on the same graph
out_long = pivot_longer(out, cols=c(2,3,4,5,6,7))
ggplot(out_long, aes(x=time, y=value, group=name, colour=name)) +
  geom_line() +
  theme_minimal() +
  xlab("Days Elapsed") +
  ylab("Number of individuals") +
  guides(colour=guide_legend(title="Compartment")) +
  ggtitle("SEIR model with homogeneity")

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
  ggtitle("SEIR Model With Population Structure: Cases")

