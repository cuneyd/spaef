function [ spaef, alpha, beta, gamma ] = SPAEF_func( sim_vec,obs_vec )

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
% [1] M. C. Demirel, J. Mai, G. Mendiguren, J. Koch, L. Samaniego, and S. Stisen, “Combining satellite data and appropriate objective functions for improved spatial pattern performance of a distributed hydrologic model,” Hydrol. Earth Syst. Sci., vol. 22, no. 2, pp. 1299–1315, Feb. 2018, doi: 10.5194/hess-22-1299-2018.
% [2] J. Koch, M. C. Demirel, and S. Stisen, “The SPAtial EFficiency metric (SPAEF): multiple-component evaluation of spatial patterns for optimization of hydrological models,” Geosci. Model Dev., vol. 11, no. 5, pp. 1873–1886, May 2018, doi: 10.5194/gmd-11-1873-2018.

  

%% correlation coefficient
cc=corrcoef(sim_vec,obs_vec);alpha=cc(1,2);

%% coefficient of variation
obs_CV=std(obs_vec)/mean(obs_vec);
sim_CV=std(sim_vec)/mean(sim_vec);
beta=sim_CV/obs_CV;

%% histogram match
obs_Zscore=zscore(obs_vec);
sim_Zscore=zscore(sim_vec);

bins=floor(sqrt(length(obs_vec)));

[N_obs,~] = histcounts(obs_Zscore,bins);
[N_sim,~] = histcounts(sim_Zscore,bins);


minOfHists = min([N_obs; N_sim], [], 1);
overlappedHist = sum(minOfHists);
histogram_match=overlappedHist/sum(N_obs);
gamma=histogram_match;

%% SPAEF
spaef = 1- sqrt( (alpha-1)^2 + (beta-1)^2 + (gamma-1)^2 );
end

