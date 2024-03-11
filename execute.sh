#!/bin/bash
set -eou pipefail

## this structure can be changed, wanted to see what information can be added
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "                           SC ANALYSIS PIPELINE             "
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

echo ""

h=`date +%H`
if [ $h -lt 12 ]; then
echo Good Morning!
elif [ $h -lt 18 ]; then
echo Good Afternoon!
else
echo Good Evening!
fi

echo ""

echo `date`

echo ""
curl -s wttr.in | head -n 7 | grep -v ┼
echo ""
echo "Authors: Louise Grimble, Dave Horsfall, Daniela Basurto-Lozada."
echo "" 
##### REMINDER TO UNCOMMENT WHEN RUNNING ON FARM #####     
#echo "Please Enter iRODS Password: "
#read -s password #-s flag used for security   
#iinit $password  

#if [ $? -ne 0 ]; then
#    echo "Password Incorrect- Please Try Again."
#else   
#    echo "Login succesful"
#fi
###################################################### 
echo ""

echo "Please Enter Path To Nextflow Script:"
read script
echo "Please Enter Nextflow Profile:"
read profile
nextflow run $script -params-file params.yaml -profile $profile -resume

echo "sc-analysis-nf pipeline completed."
echo ""