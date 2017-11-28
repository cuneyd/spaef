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

[1] Demirel, M. C., Mai, J., Mendiguren, G., Koch, J., Samaniego, L. and Stisen, S.: Combining satellite data and appropriate objective functions for improved spatial pattern performance of a distributed hydrologic model, Hydrol. Earth Syst. Sci. Discuss., 1–22, doi:10.5194/hess-2017-570, 2017a.
[2] Koch, J., Demirel, M. C. and Stisen, S.: On the importance of multiple-component evaluation of spatial patterns for optimization of earth system models &amp;amp;ndash; A case study using mHM v5.6 at catchment scale, Geosci. Model Dev. Discuss., 1–25, doi:10.5194/gmd-2017-238, 2017.
[3] Demirel, M. C., Koch, J. and Stisen, S.: SPAEF: SPAtial EFficiency, Researchgate, doi:10.13140/RG.2.2.18400.58884, 2017b.

function:
    SPAEF : spatial efficiency   
"""

# import required modules
import numpy as np
from scipy.stats import variation,zscore
######################################################################################################################
def filter_nan(s,o):
    data = np.transpose(np.array([s.flatten(),o.flatten()]))
    data = data[~np.isnan(data).any(1)]
    return data[:,0], data[:,1]
######################################################################################################################
def SPAEF(s, o, bins):
    #remove NANs    
    s,o = filter_nan(s,o)
    #compute ratio of CV
    alpha = variation(s)/variation(o)
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
    hh = np.sum(minima)/np.sum(hobs)
    #compute corr coeff
    cc = np.corrcoef(s,o)[0,1]
    #compute SPAEF finally with three vital components
    spaef = 1- np.sqrt( (cc-1)**2 + (alpha-1)**2 + (hh-1)**2 )  

    return spaef, cc, alpha, hh
######################################################################################################################
 
