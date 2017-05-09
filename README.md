# VBFSUSYSpring15

## Introduction

This analysis tool was developed to fulfill a search for Supersymmetry particles at the LHC. It is compatible with CMSSW_7_4_14 framework release and does the following things:

1) Run over 10 millions events stored in the miniaodv2 format within 10 minutes;
2) Applies 8 different event selections and creates for each one of then a collection of plots of the most important physical quantities;
3) Uses a Data-Driven method in order to predict the event distributions in the regions where statistics is low
4) Calculates cross section limits in a 3-dimensional phase-space for each different signal data samples

## Folders content

### plugins

This folder stores the core of the analysis and is divided in the following analyzers:

- VBFSUSYanalyzer: the main analyzer code. It runs different analysis selections on the same data sample and saves the histograms of the most important distributions in a .root file.

- VBFSUSYLtoTfactors: is a variation of the main analyzer. It is supposed to run exclusively on QCD data samples and applies a very loose distributions based on the properties of the reconstructed jets. The most important distributions are saved in a .root file in preparation for the data-driven background prediction techniques.

### python

 this folder stores the main .py configuration files necessary for running the analysis code though the cmsRun command (available in the CMSSW framework). Furthermore each config files stores different environment variables necessary to run different analysis selections, to run on specific data samples or enable debug mode.

### test

Stores the bash scripts necessary to interface the python config files with the job submission tool used by the NAF computing center.

### tools

Stores all the following ROOT macros:

- crvalidation: creates stacked histograms of the main background distributions of the main physical distributions for the most important control regions.

- prospects_13TeV: creates signal and background distribution predictions;

- xseclimits: is the core of the PhD thesis study. Takes all the .root files produced by the analyzers as input and searches for the cross section limit minimum versus 3 of the main analysis variables for each of the signal scenarios considered in the analysis. Furthermore creates 2-D histograms of the cross section limit distribution and makes comparison between the cross section prediction of the CMS experiment at LHC and the results coming from this analysis.
