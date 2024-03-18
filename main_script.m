%% main MATLAB script regarding the analysis of the within-host SARS-CoV-2 model 
%you can perform the whole analysis and create each figure from this script

%add path
addpath(genpath(pwd))

clearvars;
clc;
%% fitting our within-host SARS-CoV-2 model to the short-term infection dynamics by Ke et al. 2022
main_opt_Ke2022

%% plot of Figure 2 - fits and comparison to Ke et al. fits
plot_Figure2

%% plot of Figure 3 - short-term simulations
plot_Figure3

%% plot of Figure 4 - long-term fits and reinfection
plot_Figure4

%% plot of Figure 5 - reinfection simulations
plot_Figure5

%% plot of Figure 6 - reproduction number R
plot_Figure6

%% plot of Figure S1 - all model fits
plot_FigureS1

%% plot of Figure S2 - sensitivity analysis
plot_FigureS2
