# 09-25-2019 - SIR With Demographics (Birth/Death) and Stability

## SIR With Demographics

When new hosts being introduced into or removed from the population over time, we are not restricted to assuming that the pathogen eventually goes into extinction. If we want to explore dynamics of long-term persistence such as endemic or enzootic diseases, we must move past the assumption that the pathogen goes into extinction and account for demographic processes such as births (introduction of new susceptibles) and deaths. 

We can incorporate demographics into the SIR model by simply considering that hosts have a natural lifespan ($1/\mu$; which is currently at about 71 years of age). This leads to the *rate* at which hosts suffer natural mortality (does not reflect mortality by pathogen) in each compartment of the model, $\mu$. Many epidemiologists use this value to represent the birth rate as well. This is to set up conditions where the population does not grow or decline over time.

<center>

$\frac{dS}{dt} = \mu - \beta S I - \mu I$

$\frac{dI}{dt} = \beta S I - \gamma I - \mu I$

$\frac{dR}{dt} = \gamma I - \mu R$

</center>

