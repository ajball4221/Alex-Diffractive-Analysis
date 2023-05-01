#source /minerva/app/users/ajball/MINERvA101/opt/bin/setup.sh
#source /minerva/app/users/ajball/MINERvA101/opt/bin/setupROOT6OnGPVMs.sh
source /cvmfs/larsoft.opensciencegrid.org/products/setup
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
setup root v6_22_06a -q e19:p383b:prof
setup jobsub_client v1_3_3 #v_lite #v1_3_3
export JOBSUB_GROUP=minerva
source /cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/setups
setup cmake v3_9_5
source /minerva/app/users/ajball/ncdiff/opt/bin/setup.sh
source /minerva/app/users/ajball/ncdiff/CC-NuE-XSec/setup_ncdiff.sh
cd /minerva/app/users/ajball/ncdiff/
voms-proxy-destroy
kx509
voms-proxy-init --rfc --voms=fermilab:/fermilab/minerva/Role=Analysis --noregen -valid 72:7
pwd
