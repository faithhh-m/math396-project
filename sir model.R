# load required libraries
library(ggplot2)
library(deSolve)
library(tidyverse)

# plot basic SIR model using only parameters determined from dataset

# set R0 value for zika
R_0 = 3.8

# define SIR model function
sir_model <- function(times, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    N = S + I + R
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# set initial conditions, parameters, sequence of time points for solving the model
init = c(S = 7390, I = 1, R = 0)
parameters = c(beta = (R_0*0.2*7)/7391, gamma = 7/5)
times = seq(0, 60, length.out=1000)

# solve the SIR model using initial conditions, parameters and time frame defined above
out = ode(y = init, times = times, func = sir_model, parms = parameters)
out = as.data.frame(out) # convert output to dataframe for easier manipulation

# convert output from wide to long format
out_long <- pivot_longer(out, cols = c(2,3,4))
# plot the SIR model from above
ggplot(out_long, aes(x=time, y=value, group=name, colour=name)) + geom_line(lwd=2) + theme_minimal() + xlab("Time (days)") + ylab("Number of individuals") +
  guides(colour=guide_legend(title="Compartment")) + ggtitle("SIR model: beta = 7.198e-4, gamma = 0.2")
