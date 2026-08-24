seed = 20260226
jobname = bpp_run
seqfile = bpp_input_single_locus.phy
Imapfile = bpp_imap.txt
speciesdelimitation = 1 0 2
speciestree = 1 0.1 0.1 0.1
speciesmodelprior = 1
species&tree = 6 North North_Coast Sierra_1 Sierra_2 Sierra_3 South_Coast
               2 2 2 2 2 2
               (((Sierra_1:0.1333202724,Sierra_2:0.0836721975):0.2736020054,((North:0.1911468317,Sierra_3:0.1372228428):0.0583447617,North_Coast:0.1365682539):0.0402899243):1e-06,South_Coast:0.2909357973);
phase = 1 1 1 1 1 1
usedata = 1
nloci = 1
cleandata = 0
thetaprior = gamma 2 2000
tauprior = gamma 2 1000
finetune = 1
burnin = 40000
sampfreq = 10
nsample = 240000
