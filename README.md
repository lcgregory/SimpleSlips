# SimpleSlips
A simple Bayesian inference MCMC code for modelling cosmogenic data from bedrock fault scarps

This code is based on the forward model developed in Schlagenhauf et al., 2010, 'modelscarp'. Their forward model calculates the 36Cl concentration that would result from a user-defined slip history.

The MCMC code was developed by Leo Zijerveld, Patience Cowie, and Laura Gregory, first published in Cowie et al., 2017. The code published here has been modified to run faster, and with more flexibility in setting up the priors, additional scripts for calculating the posterior probability, and additional plotting scripts. The MCMC code generates potential slip histories, using a Bayesian-influened reversible jump Markov Chanin Monte Carlo approach to sample the posterior distribution. The accepted slip histories are output as an age arrary (with 'slip' as the top line of the array), with details of the liklihood and fit to the data for each model.

To run the code, download SimpleSlips and keep the same file structure on your machine. You must supply the site-specific details in the same format as required by the modelscrap code (described in Schlagenhauf et al., 2010) namely the datacolluvium.txt, datamagfield.txt, and datarock.txt. MCMC_Flex_CPSS.m should be updated with site details and a few parameters need to be set. This script is run in Matlab, and it takes ~24 hours to run 200k iterations (recommended). 

We developed various plotting scripts, which are all in SimpleSlips/plotting, and are described in the plotting readme.txt

For further information and implementation details, see:

PA Cowie, Phillips, RJ, Roberts, GP, McCaffrey, K, Zijerveld, LJJ, Gregory, LC et al. (2017). Orogen-scale uplift in the central Italian Apennines drives episodic behaviour of earthquake faults. Scientific Reports 7:44858 DOI: 10.1038/srep44858

A Schlagenhauf, Gaudemer, Y, Benedetti, L, Manighetti, I, Palumbo, L, Schimmelpfennig, I, Finkle, R, and Pou, K (2010). Uising in situ Chlorine-36 cosmonuclide to recover past earthquake histories on limestone normal fault scarps: a reappraisal of methodology and interpretations. Geophysical Journal International v 182, p 36-72, doi: 10.1111/j.1365-246X.2010.04622.x

This code was developed with funding from the UK Natural Environment Research Council (NERC) through a Standard Grant and Independent Research Fellowship.

Please contact Laura Gregory at l.c.gregory@leeds.ac.uk with questions or comments.
