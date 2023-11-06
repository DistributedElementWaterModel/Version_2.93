
#Define model domain for distributed element water balance model
stationMask control_mask.txt

#Preprocessing for distributed element water balance model
predew control_pre.txt

#Set environment variable for SeNorge meteorological data
source set_env

#Run distributed element water balance model
dew control_dew.txt

