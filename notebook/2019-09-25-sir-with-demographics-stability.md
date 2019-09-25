# 09-25-2019 - SIR With Demographics (Birth/Death) and Stability

## SIR With Demographics

When new hosts being introduced into or removed from the population over time, we are not restricted to assuming that the pathogen eventually goes into extinction. If we want to explore dynamics of long-term persistence such as endemic or enzootic diseases, we must move past the assumption that the pathogen goes into extinction and account for demographic processes such as births (introduction of new susceptibles) and deaths. 

We can incorporate demographics into the SIR model by simply considering that hosts have a natural lifespan ($1/\mu$; which is currently at about 71 years of age). This leads to the *rate* at which hosts suffer natural mortality (does not reflect mortality by pathogen) in each compartment of the model, $\mu$. Many epidemiologists use this value to represent the birth rate as well. This is to set up conditions where the population does not grow or decline over time.

<center>
  
$$\frac{dS}{dt} = \mu - \beta \text{S} \text{I] - \mu \text{I}$$

$\frac{dI}{dt} = \beta S I - \gamma I - \mu I$

$\frac{dR}{dt} = \gamma I - \mu R$

$\frac{dS}{dt} + \frac{dI}{dt} + \frac{dR}{dt} = 0$

</center>

There are some subtleties that we miss here still (newborn passive immunity, unequal rates of birth and death, etc.) but one could make the argument that these subtleties can be abstracted for the sake of practicality. There are some exceptions to this but we will move on from this for now. 

With the new set of equations that have incorporated birth/death rate, we need to consider how $R_0$ changes. $\beta$ is the transmission rate per infective and the negative terms in the equation for the change in the proportion of infectives ($-\gamma I$ & $-\mu I$) show that each infected individual spends an average of $\frac{1}{\gamma + \mu}$ time units in the infected class. If we assume that the initial population is comprised entirely of susceptibles (S=1.0), we can see that the average number of new infections per individual is equal to the transmission rate multiplied by the average infection period.

<center>

$R_0 = \frac{\beta}{\gamma + \mu}$

</center>

### The Equilibrium State

With the inclusion of the birth and death of hosts, we can start to think about persistence with a bit more detail. Are there conditions by which the host population and pathogen reach a state of equilibrium? There is an obvious one and that is the disease-free condition where the disease has gone into extinction. This condition isn't all that interesting from our point of view. This condition and another, more interesting condition can be seen if we set $\frac{dI}{dt} = 0$. 

<center>

$0 = \beta S I - \gamma I - \mu I$

$0 = I (\beta S - (\gamma + \mu))$

</center>

this statement is made true if I = 0 or S assumes a value that makes the equation zero ($\frac{\gamma + \mu}{\beta}$).


## Questions

1. Why do we account for $\mu$ twice in $\frac{dS}{dt}$?
