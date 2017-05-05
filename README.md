# VBFSUSYSpring15

## Introduction

This analysis tool is compatible with CMSSW_7_4_14 framework release and does the following things:

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

Stores all the necessary ROOT macros. Each one of the macros take one or multiple .root files produced by the analyzers (stored in the subfolders) and creates variables distribution plots, cross section limits and "signal distributions vs selection cut" predictions.
