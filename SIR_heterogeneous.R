library(ggplot2)
library(deSolve)
library(tidyverse)
library("outbreaks")

data = zika_yap_2007

# keep only the rows from day 63 and onward.
# This is because all previous weeks have 0 cases
data = data[data$nr >= 63,]

# The first case occurs on day 63, update the nr column
data$nr = data$nr - 63

# SIR model with population structure
sir_model <- function(times, state, parameters) {
    
    with(as.list(c(state, parameters)), {
        
        N = S1 + S2 + I1 + I2 + R1 + R2
        dS1 <- -beta11 * S1 * I1 - beta12 * S1 * I2
        dS2 <- -beta22 * S2 * I2 - beta21 * S2 * I1
        dI1 <-  beta11 * S1 * I1 + beta12 * S1 * I2 - gamma * I1
        dI2 <-  beta22 * S2 * I2 + beta21 * S2 * I1 - gamma * I2
        dR1 <-                                        gamma * I1
        dR2 <-                                        gamma * I2
        
        return(list(c(dS1, dS2, dI1, dI2, dR1, dR2)))
    })
}

# set initial conditions, parameters,
# sequence of time points for solving the model

init = c(S1 = 1200, S2 = 2800, I1 = 0, I2 = 3, R1 = 0, R2 = 0) 

parameters = c(beta11 = 0.5/9725,
               beta12 = 0.025/9725,
               beta21 = 0.025/9725,
               beta22 = 1/9725,
               gamma = 1/4)

times = seq(0, 150, length.out=1000)

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
    ggtitle("SIR model with homogeneity")

# We are interested in the TOTAL number of infecteds
# For each time, we should add the values for I1 and I2 in out
out$total_I = out$I1 + out$I2

out_long = pivot_longer(out, cols=c(4,5))
ggplot() +
    geom_line(data=out_long, aes(x=time, y=total_I, group=name, colour='Simulated')) +
    geom_line(data=data, aes(x=nr, y=value, color='Real Data')) +
    theme_minimal() +
    xlab("Days Elapsed") +
    ylab("Number of individuals") +
    ggtitle("SIR Model With Population Structure: Cases")

