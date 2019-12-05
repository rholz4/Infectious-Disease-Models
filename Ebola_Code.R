###################################
#                                 #
# Infectious Disease Models in R  #
#   Written By: Ryan Holzhauer    #
#                                 #
###################################

# Imagine the year is 2016. The Ebola virus outbreak in Western Africa is fading.
# The governments in this part of the world want to find a proactive solution for
# combatting a potential re-emergence of the virus. As a biological data-scientist,
# your job is to create a dynamic model that can predict the outcomes of another
# outbreak.

### Libraries----
library(deSolve)  # solves model equations
library(ggplot2)  # for creating figures
library(gridExtra)  # for making paneled figures

### Creating the Disease Model----
# With any mathematical model, it is important to start with a simplified version.
# This allows you to make sure that you defined all the core elements properly and 
# that the code functions as expected.
# The simplest form of infectious disease models include 3 classes and 2 parameters

## Classes of the Western African Population
# S <- the proportion of the population that are susceptible to the virus
# I <- the proportion of the population that are currently infected
# R <- the proportion of the population that are now removed from the population

## Model Parameters
# beta <- the probability an I individual gives the infection to an S individual
# gamma <- the rate at which an I individual dies from the infection
model_params <- c(beta = 0.05,  
                  gamma = 0.001)  # creates the parameter vector

## Initial Conditions
model_IC <- c(S = .999, 
              I = .001, 
              R = 0)  # creates the IC vector

## Time Interval
# Lets look at how the Ebola outbreak impacts the population over a the course of 
# a year
time_int <- seq(0, 365, by = 1)  # Creates sequence of numbers that create time interval

## Construct the Model
Ebola_model <- function(time, state, pars) { 
  with(as.list(c(state, pars)), {   # defines the model as a solveable function
    dS <- 0 - beta * S * I  # describes the susceptible/"healthy" class
    dI <- beta * S * I - gamma * I  # describes the infected class
    dR <- gamma * I # describes the removed class
    return(list(c(dS, dI, dR)))  # creates a solution vector  
  })
}

## Solve the Model
model_solution <- as.data.frame(ode(func = Ebola_model,  
                                    y = model_IC, 
                                    parms = model_params, 
                                    times = time_int))
# organizes the solutions found by the ode package into the form of a data frame.
# this makes graphing the model outputs much easier

### Visualizing the Model----
# Creating a plot
(Ebola_plot <- ggplot(model_solution) +  # creates a plot using our solution data set
                     geom_line(aes(x = time, y = S, color = "S")) +  
                     geom_line(aes(x = time, y = I, color = "I")) +
                     geom_line(aes(x = time, y = R, color = "R")) +  # makes line for
                     # each class of the model and adds proper label to legend
                     scale_x_continuous(limits = c(0,370),
                                        breaks = seq(0,370,50),  # adjusts scale of
                    # x axis
                                        expand = c(0,0)) +  # removes unnecessary space
                     scale_y_continuous(limits = c(0,1), 
                                        breaks = seq(0,1,.2), 
                                        expand = c(0,0)) +
                     theme_classic() +  # creates a clear background
                     labs(color="Legend",  # adds color-coded legend
                         x="\nTime (Days)",  # creates x-axis label
                         # \n adds space in between x-axis and the label
                         y="Proportion of Population\n",  # creates y-axis label
                         title = "Impact of Ebola Outbreak on West Africa"))  # creates title
# Saving the plot
ggsave(Ebola_plot, filename = "Ebola_Plot.png",  # saves plot in this file type
       height = 5, width = 7,  # sets dimensions
       path = "Images")  # saves to proper folder

### Creating an improved Model----
# We make a lot of assumptions with our previous model that are not biologically
# accurate. One way to improve a model is to try to account for some of these assumptions.
# It is important to note that you can never account for EVERY factor within a system.

## Classes of the Western African Population
# One way to improve the relevancy of our model is to add more classes. In the previous
# model we assume that a healthy individual immediately becomes infected when they come
# in contact with the virus. Generally that is not the case in disease models. We can
# introduce an "exposed" class. In additon the current model assumes that once you
# contract Ebola you will eventually die. We can account for this by adding a "quarantine"
# class

# S <- the proportion of the population that are susceptible to the virus
# E <- the proportion of the population that have contracted the virus but are not
# yet contagious
# I <- the proportion of the population that are currently infected
# Q <- the proportion of the population that contracted the virus and have been placed
# into some form of medical care
# R <- the proportion of the population that are now removed from the population

## Model Parameters
# Since we added new classes we need to add parameters to establish how they interact
# within our model

# beta <- the probability an I individual contracts the infection to an S individual
# delta <- the rate at which an E individual becomes an I individual
# gamma <- the rate at which an I individual dies from the infection
# alpha <- the percentage of I individuals taken into medical care
# kappa <- the rate at which a Q individual becomes an S individual

improved_params <- c(beta = 0.05,
                  delta = 1/14,
                  gamma = 0.001,
                  alpha = .0001,
                  kappa = 1/30)  # creates the parameter vector

## Initial Conditions
# We will just assume the E and Q classes take up 0% of the population at the start
improved_IC <- c(S = .999, 
                 E = 0, 
                 I = .001, 
                 Q = 0, 
                 R = 0)  # creates the IC vector

improved_Ebola_model <- function(time, state, pars) { 
  with(as.list(c(state, pars)), {   # defines the model as a solveable function
    dS <- 0 - beta * S * I + kappa * Q  # susceptible class
    dE <- beta * S * I - delta * E  # exposed class
    dI <- delta * E - gamma * I - alpha * I  # infected class
    dQ <- alpha * I - kappa * Q  # quarantine class
    dR <- gamma * I  # removed class
    return(list(c(dS, dE, dI, dQ, dR)))  # creates a solution vector  
  })
}

## Solve the Model
improved_solution <- as.data.frame(ode(func = improved_Ebola_model,  
                                    y = improved_IC, 
                                    parms = improved_params, 
                                    times = time_int))
# organizes the solutions found by the ode package into the form of a data frame.
# This makes graphing the model outputs much easier

### Visualizing the Improved Model----
# Creating a plot
(Improved_plot <- ggplot(improved_solution) +  
   geom_line(aes(x = time, y = S, color = "S")) +
   geom_line(aes(x = time, y = E, color = "E")) +
   geom_line(aes(x = time, y = I, color = "I")) +
   geom_line(aes(x = time, y = Q, color = "Q")) +
   geom_line(aes(x = time, y = R, color = "R")) +  
   scale_x_continuous(limits = c(0,370),
                      breaks = seq(0,370,50), 
                      expand = c(0,0)) +
   scale_y_continuous(limits = c(0,1), 
                      breaks = seq(0,1,.2), 
                      expand = c(0,0)) +
   theme_classic() +  
   labs(color="Legend",  
        x="\nTime (Days)",  
        y="Proportion of Population\n",  
        title = "Impact of Ebola Outbreak on West Africa"))

# Saving the plot
ggsave(Improved_plot, filename = "Improved_Ebola_Plot.png",
       height = 5, width = 7, 
       path = "Images")  

### Comparing our Two Models
# One possible way to compare the results of the two models is to create a paneled
# figure

# Combining the two plots
Ebola_model_comparison <- grid.arrange(Ebola_plot, Improved_plot, ncol = 2)
# creates a paneled figure with 2 columns

# Saving the figure
model_graphic <- arrangeGrob(Ebola_plot, Improved_plot, ncol = 2)  # stores the 
# paneled figure as an element that can be saved
ggsave(model_graphic, filename = "Ebola_Model_Comparison.png",
       height = 4,width = 9, 
       path = "Images")  
