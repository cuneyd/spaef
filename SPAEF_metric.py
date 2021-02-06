#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:33:33 2017
@ authors:                 Mehmet Cüneyd Demirel, Gorka Mendiguren, Julian Koch, Simon Stisen and Fabio Oriani
@ author's website:        http://www.space.geus.dk/
@ author's webpage:        http://akademi.itu.edu.tr/demirelmc/
@ author's email id:       demirelmc@itu.edu.tr

A libray with Python functions for calculation of spatial efficiency (SPAEF) metric.

Literature:

[1] Demirel, M. C., Mai, J., Mendiguren, G., Koch, J., Samaniego, L., & Stisen, S. (2018). Combining satellite data and appropriate objective functions for improved spatial pattern performance of a distributed hydrologic model. Hydrology and Earth System Sciences, 22(2), 1299-1315. https://doi.org/10.5194/hess-22-1299-2018
[2] Koch, J., Demirel, M. C., & Stisen, S. (2018). The SPAtial EFficiency metric (SPAEF): multiple-component evaluation of spatial patterns for optimization of hydrological models. Geoscientific Model Development, 11(5), 1873-1886. https://doi.org/10.5194/gmd-11-1873-2018

Cite as: Demirel, M.C., Koch, J., Stisen, S., 2018. SPAEF: SPAtial EFficiency. GitHub. https://doi.org/10.5281/ZENODO.1158890

function:
    SPAEF : spatial efficiency   
"""

# import required modules
import numpy as np
from scipy.stats import variation,zscore
import math
######################################################################################################################
def filter_nan(s,o):
    data = np.transpose(np.array([s.flatten(),o.flatten()]))
    data = data[~np.isnan(data).any(1)]
    return data[:,0], data[:,1]
######################################################################################################################
def SPAEF(s, o):
    #remove NANs    
    s,o = filter_nan(s,o)
    
    bins=np.around(math.sqrt(len(o)),0)
    #compute corr coeff
    alpha = np.corrcoef(s,o)[0,1]
    #compute ratio of CV
    beta = variation(s)/variation(o)
    #compute zscore mean=0, std=1
    o=zscore(o)
    s=zscore(s)
    #compute histograms
    hobs,binobs = np.histogram(o,bins)
    hsim,binsim = np.histogram(s,bins)
    #convert int to float, critical conversion for the result
    hobs=np.float64(hobs)
    hsim=np.float64(hsim)
    #find the overlapping of two histogram      
    minima = np.minimum(hsim, hobs)
    #compute the fraction of intersection area to the observed histogram area, hist intersection/overlap index   
    gamma = np.sum(minima)/np.sum(hobs)
    #compute SPAEF finally with three vital components
    spaef = 1- np.sqrt( (alpha-1)**2 + (beta-1)**2 + (gamma-1)**2 )  

    return spaef, alpha, beta, gamma
######################################################################################################################
 
