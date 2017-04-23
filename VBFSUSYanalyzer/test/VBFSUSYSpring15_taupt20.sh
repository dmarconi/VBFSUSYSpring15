#!/bin/sh

source /afs/desy.de/user/d/dmarconi/VBF-LS-tau/setenv.sh

cmsRun /nfs/dust/cms/user/dmarconi/workdir/VBFsignalSamples_Madgraph/CMSSW_7_4_14/src/VBFSUSYSpring15/VBFSUSYanalyzer/python/ConfFile_gridcontrol_taupt20_cfg.py $1 $2
