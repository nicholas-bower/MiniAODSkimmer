# MiniAODSkimmerThis package is used to apply cleaning on the MiniAOD samples. We do the Electron/Muon Cleaning in this step and store results for Taus  collections cleaned. # Setting up the environment:```bash$ setenv SCRAM_ARCH slc7_amd64_gcc900 $ cmsrel CMSSW_12_0_0_pre4$ cd CMSSW_12_0_0_pre4/src$ git cms-init$ cmsenv$ git clone --recursive git@github.com:red1habibullah/MiniAODSkimmer.git -b UL_12X_2018$ git cms-addpkg PhysicsTools/PatAlgos$ rm PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc $ cp /uscms_data/d3/rhabibul/MiniAODCleaner/CMSSW_10_6_20/src/PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc$ scram b -j8```