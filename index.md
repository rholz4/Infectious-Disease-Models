---
layout: post
title: "Introduction to Infectious Disease Models"
date: 2019-12-5 11:02:30
author: Ryan Holzhauer
meta: "Tutorials"
---

### Tutorial Aims

#### <a href="#section1"> 1. Understand the concept of an infectious disease model</a>

#### <a href="#section2"> 2. Construct and solve a model using deSolve</a>

#### <a href="#section3"> 3. Visualize and interpret model results</a>

#### <a href="#section4"> 4. Challenge: Explore ways to improve models</a>

<center><img src="{{ site.baseurl }}/Infectious_Bacteria.jpg" alt="Img"></center>

In this tutorial we are going to use infectious disease models to better understand
the dynamics behind a potential Ebola virus outbreak. Disease models are a powerful
predictive tool that allows scientists and mathematicians to forecast the effects
of an infection within a given population. These models predict how much of the
population will be affected by an infection, how fast the infection will spread
and are used to optimize potential treatment options.
---------------------------

You can get all of the resources for this tutorial from <a href="https://github.com/rholz4/Infectious_Disease_Models" target="_blank">this GitHub repository</a>. Clone and download the repo as a zip file, then unzip it._

<a name="section1"></a>

## 1. Understand the concept of an infectious disease model

### What is mathematical modelling?
In simplest terms, a mathematical model is the translation of a real-life problem
into a series of interconnected mathematical formulas. Often times in the real world
many problems can be incredibly difficult to solve solely using logic or information
within a given field. From analyzing how introducing a new species will effect its
ecosystem or when to invest into a particular stock, mathematical models provide a
powerful computational tool for analyzing these problems.

### What is an Infectious Disease Model?
Many mathematical models that are used to describe biological and ecological problems
are often classified as dynamic. This means that the elements within the system
function in a way that is dependent on time. In mathematical terms, this dynamic
relationship can be accounted for by including derivatives in the equations.
The proper name for an equation that contains a derivative is a _differential equation_.
An infectious disease model is a type of dynamic model that describes how a certain
disease spreads throughout a population. These models are often organized by grouping
individuals into categories based on how they are impacted by the disease (ex. healthy,
infected, dead, etc.). These groupings are referred to as _classes_.

### How do These Models Work?
In concept dynamic models function very similarly to Bayesian statistical models.
Both types of models utilize prior information to extrapolate results for their given
systems. In dynamic models these elements are referred to as _parameters_

There are a couple key differences though. Instead of using a data set as the element
you are trying to analyze, dynamic models use mathematical equations. Dynamic models
are solved by taking elements referred to as _initial conditions_, inserting them
into the corresponding equations, solving the equations individually, and repeating this
process over a determined sequence of time. In other words the model is solved algorithmically.
This is different to the iterative approach used to solve Bayesian models.

### Why do We Need Coding/Programming Analyze Infectious Disease Models?
As you would expect solving a dynamic model by hand is not only a time-consuming
process but also very difficult. Technological advancements over the last 20 years
have made these problems much easier to analyze and refine. As a result, this niche
area of biomathematical research is rapidly growing.

Multiple different programming languages can be used to solve Infectious Disease
Models. For the sake of this tutorial we will use R since it also allows us to create
visually-appealing figures that display the model results.

<a name="section2"></a>

## 2. Construct and Solve an Infectious Disease Model
Before we start conceptualizing the model, open up `RStudio` and create a new R script.
Include your name, the date and the product you are trying to produce at the top
of the script. You can use the "#" key to convert the title of the script into a comment.
This way the text will not interfere with the actual code. Once this is done you
can then load the packages needed to perform this tutorial
```r
# Libraries
library(deSolve)  # solves model equations
library(ggplot2)  # for creating figures
library(gridExtra)  # for making paneled figures
```
You can also set the working directory of the script if you want to work within a
specific folder of your repository.

### Defining the Scenario and the Research Question
Imagine the year is 2016. The Ebola virus outbreak in Western Africa is fading.
The governments in this part of the world want to find a proactive solution for
combatting a theoretical re-emergence of the virus. As a talented independent
researcher, you are tasked to create an infectious disease model that can be used
to predict the outcomes of another outbreak.

We can define the research question to be "How will a re-emergence of the Ebola virus
affect the population of Western Africa?". While this example seems a bit hyperbolic,
this is relatively how similar projects are started by scientists today. Generally
after an infection or virus reaches a population, scientists find out via news outlets.
If they are interested/motivated by the problem, they then go about researching and
studying it.

### Creating the Infectious Disease Model
With any mathematical model, it is important to start with a simplified version.
This allows you to make sure that you defined all the core elements properly and
that the code functions as expected. The simplest form of infectious disease models
include 3 classes and 2 parameters and are referred to as _SIR models_

#### Classes of the Western African Population
The first step of making an infectious disease model is to define the different
classes within the model. In the most simple model populations are grouped into
three different classes: susceptible, infected and removed (this is where the SIR
abbreviation comes from).
```r
### Constructing the Model----
## Defining the classes
# S <- the proportion of the population that are susceptible to the virus
# I <- the proportion of the population that are currently infected
# R <- the proportion of the population that are now removed from the population
```
These classes do not actually need to be defined in the code but it is helpful to
include them in the script so that way it is easier to interpret results later on.

#### Model Parameters
After the classes are established, we need to define how the Ebola virus travels
through the three classes. This can be explained within the model by parameters.
In dynamic models a parameter can either represent the probability that an event
takes place or the rate that an event takes place.

We know that the Ebola virus spreads when an infected individual comes into contact
with a healthy individual. In addition we know that someone infected with Ebola can
live for a certain length of time before they die. We can define both of these ideas
using parameters
```r
# beta <- the probability an I individual gives the infection to an S individual
# gamma <- the rate at which an I individual dies from the infection.
model_params <- c(beta = 0.05,
                  gamma = 0.001)  # creates the parameter vector
```
We will assume for now that the values for these parameters are somewhat accurate.
The parameters were stored within a vector so that they can be easily accessed when
we are ready to solve the model.

#### Initial Conditions
The other key piece of prior information we need is the initial conditions of the
model. In our case we want our initial conditions to represent the proportion of the
Western African population that falls within each class at the start of the outbreak.
```r
model_IC <- c(S = .999,
              I = .001,
              R = 0)  # creates the IC vector
```
Since the population of Western Africa is around 350,000,000 people, these conditions
implied that we are starting our model after 350,000 people had already been infected.

#### Time Interval
From our research question we want to look at how the Ebola outbreak impacts the
population over a the course of a year in time. Instead of creating a vector to
represent the time interval, we can just create a numerical sequence.
```r
time_int <- seq(0, 365, by = 1)  # Creates sequence of numbers that create time interval
```
In this case this sequence produces contains every whole number from 0 to 365.

#### Construct the Model
Now that we have all the necessary components, we can write the model. There are multiple
ways to write and solve a dynamic model within R but we will choose to focus on the
one enabled by the `deSolve` package. For this approach we just need to define the
model as a function of its three mathematical components: time, state and parameters.

A great way to to start writing a model is to create a something called a _wire diagram_.
These diagrams are helpful for visualizing how the parameters interact with each class.

Insert wire diagram image here

In a wire diagram, each arrow represents an event that occurs in the system and the
direction in which it travels. If the arrow enters a class, then we add the corresponding
expression in its equation. If the arrow leaves a class, then we subtract.

Now that we understand how the events travel through the syste, we just need to
create the mathematical expressions for them. The easiest way to do this is to look
back at how we defined the parameters.

We defined _beta_ as the __probability__ that the infection is spread when an I
individual __comes into contact__ with an S individual. The action of contact implies
that the S and I populations have to interact in order to spread the virus. Similarly to
in physics, this can be described by S * I. Since we know that this event is not guaranteed
to happen, we can further describe the transmission of the virus as beta * S * I.
We defined _gamma_ as the __rate__ at which I individuals enter the R class. This can
simply be described by gamma * I.

We can now write the function that represents the model.
```r
Ebola_model <- function(time, state, pars) {
  with(as.list(c(state, pars)), {   # defines the model as a solveable function
    dS <- 0 - beta * S * I  # describes the susceptible/"healthy" class
    dI <- beta * S * I - gamma * I  # describes the infected class
    dR <- gamma * I # describes the removed class
    return(list(c(dS, dI, dR)))  # creates a solution vector
  })
}
```
Here we define each equation by adding a "d" to the front of each class. This indicates
that these are differential equations.
The 'return' argument at the end of the function is what allows the model to update
its solution after every step in the time interval.

#### Solve the Model
The last step is to solve the model. We will use the `ode` function within the `deSolve`
package to solve all the differential equations within the model. This function's
outputs can then be organized into a data frame so that we can construct figrues after.
```r
model_solution <- as.data.frame(ode(func = Ebola_model,
                                    y = model_IC,
                                    parms = model_params,
                                    times = time_int))
# organizes the solutions found by the ode package into the form of a data frame.
# this makes graphing the model outputs much easier
```

## 3. Visualize and Interpret Model Results

### Creating a Line Plot for our Model
We want to show how the classes of the Western African population changed over the
course of the Ebola outbreak. The best way to visualize this is to generate a line
plot using the `The ggplot2` package. This allows us to manipulate many of the different
elements of the figure.
For more on `ggplot2`, read the official <a href="https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf" target="_blank">ggplot2 cheatsheet</a>.
```r
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
```
The resulting graph should look like this:
<center> <img src="{{ site.baseurl }}/Ebola_Plot.png" alt="Img" style="width: 800px;"/> </center>

### Interpreting the Line Plot
When you are trying to interpret the results of an infectious disease model, it is
__really__ important to keep in mind the prior information established when writing
the model. It is obvious from this graph and observe that after about 9 months the
West African population will all be infected with Ebola or dead. What you really
want to be looking to address is __why__ this result is observed.

One potential explanation
to this question would be to look at the parameters. The beta value (0.05) is significantly
larger than the gamma value (0.001). This means that the virus defined in the model
is not only spreading at a relatively high probability (a value of 0.02 is generally
considered to be at an "epidemic" level) but the I individuals remain in the population
for a very long time. The increase of I individuals in turn makes it more likely for
an S individual to come into contact and contract Ebola from an I individual. Because
of this occurrence, logistic decay/growth is observed in the S and I class respectively.

Another potential explanation could be that the initial proportion of I individuals is
too large in our model. If you consider that the population of West Africa is about
350 million people then the model suggests that there are 350,000 infected individuals
at day 0. This is well past the introduction of "patient zero" into the population.
If the initial proportions were shifted to be more biologically accurate it is likely
that the outbreak would be much less detrimental over the course of the first year.

At this point we have now gone over how to go about creating an infectious disease
model. In contrast to statistical modelling, the analysis and interpretations of the
results are easier to understand and explain. From my experience the most challenging
part of making an infectious disease model is determining the prior information.
Often times you have very little data or background knowledge about the system you
are trying to model. This makes finding the numeric value for certain parameters (ex.
how infectious a disease is) very difficult. Parameters can be very sensitive; meaning
that a small change in its value could drastically alter the behavior of the model.
To try and mitigate this problem many recent studies try to use an machine learning
approach to estimate parameters. More information about this idea can be found in
__this article__


<a name="section4"></a>

## 4. Bonus: Explore Ways to Improve Models
From our interpretations of the first model we point out some of the different flaws
that make the results unrealistic. There are _many_ different ways to approach improving
our simple model. It is important to keep in mind though that __no model is perfect__;
meaning that you can never create an exact replica of a biological event. There are
countless factors that could contribute to how an infectious disease model functions
that and it is impossible to account for every single one _and_ have the model still
function properly. You also want to make sure that the refined model is still __biologically
relevant__. Often times researchers that develop mathematical models will instinctively
alter the model so that it becomes easier to analyze. It is important to always make
sure you can justify changes to the model based on what is you know about the biology
of the problem.

In this section I will walk through an example of how I would go about making the
Ebola model more realistic. Note that this example is by no means optimized and is
meant to just showcase the logic behind some of the changes.  

### Population Classes
One way to improve the relevancy of this model would be to add more classes. In the
previous model we only separate the population of Western Africa into a healthy,
infected and removed classes. One flaw in this framework is that this class structure
assumes that a healthy individual immediately becomes infected when they come in
contact with the virus. Generally that is not the case with infectious diseases.
To account for this we can introduce an "exposed" class (E). In additon, the current
model assumes that once you contract Ebola you will eventually die and there is no
way to treat it. We can account for this by adding a "quarantine" class (Q).

```r
### Creating an improved Model----
## Classes
# S <- the proportion of the population that are susceptible to the virus
# E <- the proportion of the population that have contracted the virus but are not
# yet contagious
# I <- the proportion of the population that are currently infected
# Q <- the proportion of the population that contracted the virus and have been placed
# into some form of medical care
# R <- the proportion of the population that are now removed from the population
```

### Parameters
Since we have added two new classes to the model we need to add at least two more
parameters that interconnect these classes. Here is one potential way to do this:

```r
## Model Parameters
# beta <- the probability an I individual transmits the infection to an S individual
# delta <- the rate at which an E individual becomes an I individual
# gamma <- the rate at which an I individual dies from the infection
# alpha <- the percentage of I individuals taken into medical care
# kappa <- the rate at which a Q individual becomes an S individual

improved_params <- c(beta = 0.05,
                  delta = 1/14,
                  gamma = 0.001,
                  alpha = .0001,
                  kappa = 1/30)  # creates the parameter vector
```

### Initial Conditions
We will assume that the population has the same initial conditions as the previous
model so that way we can directly compare how adding the E and Q classes affect the
behavior of the model.
```r
## Initial Conditions
# We will assume the E and Q classes take up 0% of the population at the start
improved_IC <- c(S = .999,
                 E = 0,
                 I = .001,
                 Q = 0,
                 R = 0)  # creates the IC vector
```

### Constructing the model
Using the prior information established above, try and write the code for the model
yourself. A good first step would be to draw a wire diagram to help visualize what
is going on in the model. It is probably easiest to draw it on a piece of paper;
that way you can easily make adjustments. Once you have drawn the wire diagram, modify
the code used for the simple model by adding in the newly-defined prior information.
Once you have constructed the model, try and solve it.

_Hint 1: When drawing the wire diagram, think about how each parameter is defined.
What classes does each parameter affect? How do they affect each class?_
_Hint 2: Think about how the simple model was constructed. Which parameters are rates
and which are probabilities? How do we know if we are adding or subtracting expressions
from an equation?_

Solution:
```r
## Constructing the model
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
```

### Visualizing the Model
Once you have created and solved the improved model, you can now create a line plot.

```r
### Visualizing the Improved Model----
## Creating a plot
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

## Saving the plot
ggsave(Improved_plot, filename = "Improved_Ebola_Plot.png",
       height = 5, width = 7,
       path = "Images")
```
The resulting graph should look like this:
<center> <img src="{{ site.baseurl }}/Improved_Ebola_Plot.png" alt="Img" style="width: 800px;"/> </center>

### Comparing our Two Models
Now that we have plots for both of the Ebola models, we can create a paneled figure
so we can easily compare them. To do this we will utilize the package `gridExtra`.
```r
### Comparing the Simple and Improved Ebola Models----
## Combining the two plots
Ebola_model_comparison <- grid.arrange(Ebola_plot, Improved_plot, ncol = 2)
# creates a paneled figure with 2 columns

## Saving the figure
model_graphic <- arrangeGrob(Ebola_plot, Improved_plot, ncol = 2)  # stores the
# paneled figure as an element that can be saved
ggsave(model_graphic, filename = "Ebola_Model_Comparison.png",
       height = 4,width = 9,
       path = "Images")
```
The resulting graph should look like this:
<center> <img src="{{ site.baseurl }}/Ebola_Model_Comparison.png" alt="Img" style="width: 800px;"/> </center>

Using the paneled figure, how does the behavior of the models differ? Do both models
predict the same outcome for the outbreak? Can you extrapolate why the classes in
the improved model behave the way they do? What can you do to improve the model
even more?

## Conclusion
_In this tutorial we learned:_
#### 1. What an infectious disease model is

#### 2. How to construct and solve an infectious disease model

#### 3. How to visualize and interpret dynamic model results

#### 4. Different ways to improve models

<hr>
<hr>

#### Check out our <a href="https://ourcodingclub.github.io/links/" target="_blank">Useful links</a> page where you can find loads of guides and cheatsheets.

#### If you have any questions about completing this tutorial, please contact us on ourcodingclub@gmail.com

#### <a href="INSERT_SURVEY_LINK" target="_blank">We would love to hear your feedback on the tutorial, whether you did it in the classroom or online!</a>

<ul class="social-icons">
	<li>
		<h3>
			<a href="https://twitter.com/our_codingclub" target="_blank">&nbsp;Follow our coding adventures on Twitter! <i class="fa fa-twitter"></i></a>
		</h3>
	</li>
</ul>

### &nbsp;&nbsp;Subscribe to our mailing list:
<div class="container">
	<div class="block">
        <!-- subscribe form start -->
		<div class="form-group">
			<form action="https://getsimpleform.com/messages?form_api_token=de1ba2f2f947822946fb6e835437ec78" method="post">
			<div class="form-group">
				<input type='text' class="form-control" name='Email' placeholder="Email" required/>
			</div>
			<div>
                        	<button class="btn btn-default" type='submit'>Subscribe</button>
                    	</div>
                	</form>
		</div>
	</div>
</div>
