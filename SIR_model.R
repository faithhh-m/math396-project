library("outbreaks")
library(ggplot2)
library(deSolve)
library(tidyverse)

data = zika_yap_2007

# keep only the rows from day 63 and onward.
# This is because all previous weeks have 0 cases
data = data[data$nr >= 63,]

# The first case occurs on day 63, update the nr column
data$nr = data$nr - 63

# Plot the cases per week
ggplot(data, aes(x=nr, y=value)) +
    geom_point(size=2, color='red') +
    geom_line() +
    xlab("Days Elapsed") + 
    ylab("Cases") +
    ggtitle("Cases Per Week During 2007 Zika Outbreak in Yap")

# basic SIR model using parameters determined from dataset
sir_model <- function(times, state, parameters) {
    
    with(as.list(c(state, parameters)), {
        
        N = S + I + R
        dS <- -beta * S * I
        dI <-  beta * S * I - gamma * I
        dR <-                gamma * I
        
        return(list(c(dS, dI, dR)))
    })
}

init = c(S = 7390, I = 1, R = 0)
parameters = c(beta = 1/9725, gamma = 1/5)
times = seq(0, 70, length.out = 1000)

out = ode(y = init, times = times, func = sir_model, parms = parameters)
out = as.data.frame(out)

out_long <- pivot_longer(out, cols = c(2,3,4))
ggplot(out_long, aes(x=time, y=value, group = name, colour = name)) +
    geom_line() +
    theme_minimal() +
    xlab("Time (days)") +
    ylab("Number of individuals") +
    guides(colour=guide_legend(title="Compartment")) +
    ggtitle("SIR model")

# plot I curve using the data above
# Did not end up including this plot in the report
infected = pivot_longer(out, cols = c(3))
ggplot() +
    geom_line(data=infected, aes(x=time, y=value, group=1, color='SIR model')) +
    geom_line(data=data, aes(x=nr, y=value, group=1, color='Real Data')) +
    theme_minimal() +
    xlab("Time") +
    ylab("Number of infected individuals") +
    guides(colour=guide_legend(title="Compartment")) +
    ggtitle("SIR model vs Real Data")
