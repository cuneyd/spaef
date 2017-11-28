function [ spaef ] = SPAEF_func( sim_vec,obs_vec,bins )

% This is a script to calculate spatial efficiency (SPAEF) which can be used to compare spatial patterns in two raster maps 
% Note that the NaN values in 2D maps are cleared first. Observed and Simulated two vectors (1D) are provided to the function.

% Detailed explanation goes here

% The newly proposed spatial efficiency metric (SPAEF) is proven to be robust 
% and easy to interpret due to its three distinct and complementary components of correlation, variance and histogram matching.

% Created on Thu Sep 13 11:33:33 2017
% @ authors:                 Mehmet Cüneyd Demirel, Gorka Mendiguren, Julian Koch, Simon Stisen and Fabio Oriani
% @ author's website:        http://www.space.geus.dk/
% @ author's webpage:        http://akademi.itu.edu.tr/demirelmc/
% @ author's email id:       demirelmc@itu.edu.tr
% 
% A libray with Matlab functions for calculation of spatial efficiency (SPAEF) metric.
% 
% Literature:
% 
% [1] Demirel, M. C., Mai, J., Mendiguren, G., Koch, J., Samaniego, L., & Stisen, S. (2017). Combining satellite data and appropriate objective functions for improved spatial pattern performance of a distributed hydrologic model. Hydrology and Earth System Sciences Discussions, 1–22. https://doi.org/10.5194/hess-2017-570
% [2] Koch, J., Demirel, M. C., & Stisen, S. (2017). On the importance of multiple-component evaluation of spatial patterns for optimization of earth system models - A case study using mHM v5.6 at catchment scale. Geoscientific Model Development Discussions, 1–25. https://doi.org/10.5194/gmd-2017-238
% 
  

%% correlation coefficient
cc=corrcoef(sim_vec,obs_vec);cc=cc(1,2);

%% coefficient of variation
obs_CV=std(obs_vec)/mean(obs_vec);
sim_CV=std(sim_vec)/mean(sim_vec);
alpha=sim_CV/obs_CV;

%% histogram match
obs_Zscore=zscore(obs_vec);
sim_Zscore=zscore(sim_vec);

[N_obs,~] = histcounts(obs_Zscore,bins);
[N_sim,~] = histcounts(sim_Zscore,bins);


minOfHists = min([N_obs; N_sim], [], 1);
overlappedHist = sum(minOfHists);
histogram_match=overlappedHist/sum(N_obs);

%% SPAEF

spaef = 1- sqrt( (cc-1)^2 + (alpha-1)^2 + (histogram_match-1)^2 );
end

