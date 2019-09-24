# 09-24-2019 - Equilibrium Stability Analysis and Next Generation Method

## Model Dynamics

As parameters vary, a model may predict different outcomes. Flexibility in our approach may be needed. One may ask how we can anticipate trajectories without resorting to extensive numerical integration? 

Under what conditions will an infectious disease invade a system?

## Invasion Threshold

What are the initial conditions that will determine whether the infection will be an epidemic or fail to spread through the population? We can answer this question by rewriting the equation for the change in proportion of infecteds over time. 

A successful invasion will only take place if the change in the proportion of infecteds is greater than zero. This change in infecteds over time might be less than zero if the initial fraction of susceptibles (S(0)) is less than the recovery rate ($\gamma$) divided by the infection rate ($\beta$). This is what we refer to as the "threshold phenomena." This is because the proportion of susceptibles in the population must exceed this critical threshold for an infection to invade. We could also think about it in relation to the relative removal rate ($\gamma/\beta$). If the relative removal rate is too large, the infection will fade away. However, this relative removal rate can be small enough to allow the infection to spread without being removed too quickly.

The relative removal rate leads us to one of the most important concepts or quantities in epidemiology. The inverse of the relative removal rate ($\beta/\gamma$) is the __basic reproduction ratio__ (R<sub>0</sub>) which is defined as the average number of secondary cases arising from an average primary case in an entirely susceptible population. In a closed population, an infectious disease with a specified R<sub>0</sub> can invade only if there is a threshold fraction of susceptibles greater than $1/R_0$