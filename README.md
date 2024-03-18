# COVID_model
contains code accompanying the paper: 
## A mathematical model for the within-host (re)infection dynamics of SARS-CoV-2
Lea Schuh<sup>a,\*</sup>, Peter V. Markov<sup>a,b</sup>, Vladimir M. Veliov<sup>c</sup>, Nikolaos I. Stilianakis<sup>a,d,\*</sup> \
<sup>a</sup> Joint Research Centre (JRC), European Commission, Via Enrico Fermi 2749, Ispra, 21027, Italy \
<sup>b</sup> London School of Hygiene & Tropical Medicine, University of London, Keppel Street, London, WC1E 7HT, United Kingdom \
<sup>c</sup> Institute of Statistics and Mathematical Methods in Economics, Vienna University of Technology, Wiedner Hauptstraße 8-10, Vienna, 1040, Austria  
<sup>d</sup> Department of Biometry and Epidemiology, University of Erlangen–Nuremberg, Waldstraße 6, Erlangen, 91054, Germany \
<sup>\*</sup> corresponding authors \
\
e-mails: lea.schuh@ec.europa.eu, nikolaos.stilianakis@ec.europa.eu

### Summary
Interactions between SARS-CoV-2 and the immune system during infection are complex. However, understanding the within-host SARS-CoV-2 dynamics is of enormous importance for clinical and public health outcomes. Current mathematical models focus on describing the within-host SARS-CoV-2 dynamics during the acute infection phase. Thereby they ignore important long-term post-acute infection effects. We present a mathematical model, which not only describes the SARS-CoV-2 infection dynamics during the acute infection phase, but extends current approaches by also recapitulating clinically observed long-term post-acute infection effects, such as the recovery of the number of susceptible epithelial cells to an initial pre-infection homeostatic level, a permanent and full clearance of the infection within the individual, immune waning, and the formation of long-term immune capacity levels after infection. Finally, we used our model and its description of the long-term post-acute infection dynamics to explore reinfection scenarios differentiating between distinct variant-specific properties of the reinfecting virus. Together, the model’s ability to describe not only the acute but also the long-term post-acute infection dynamics provides a more realistic description of key outcomes and allows for its application in clinical and public health scenarios. \
\
Required software: MATLAB (R2023a) 

### Folder structure
<ul>
  <li>Data</li>
  <li>Figures</li>
  <li>Results</li>
  <li>Scripts</li>
</ul> 

### How to run the code
The important script is main_script.m for running the workflow (parameter estimation, simulation of model fits of Ke et al. 2022, comparisons, simulation study, long-term infection simulations, and reinfection simulations) and for plotting the resulting figures. 

#### Data
Nasal CN values of 60 individuals as measured by Ke et al., 2022, and the individual-specific estimated parameters of Ke et a., 2022, used for comparison. 

#### Figures
All raw and final subfigures and figures in png, pdf, and svg formats. 

#### Results
The individual-specific estimated parameters using our model. 

#### Scripts
All functions used by the main_script.m can be found here. More information on the functions can be found in the scripts themselves eg. input, output.
