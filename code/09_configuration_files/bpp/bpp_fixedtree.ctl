seed = 20260226
jobname = bpp_fixedtree
seqfile = bpp_input_single_locus.phy
Imapfile = bpp_imap.txt
speciesdelimitation = 0
speciestree = 0
speciesmodelprior = 1
species&tree = 6 North North_Coast Sierra_1 Sierra_2 Sierra_3 South_Coast
               2 2 2 2 2 2
               ((North,North_Coast),((Sierra_1,Sierra_2),(Sierra_3,South_Coast)));
phase = 1 1 1 1 1 1
usedata = 1
nloci = 1
cleandata = 0
thetaprior = gamma 2 2000
tauprior = gamma 2 1000
finetune = 1
print = 1 0 0 0
burnin = 5000
sampfreq = 10
nsample = 50000
