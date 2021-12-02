#! /bin/bash
#this is run simulations on pc3 and store output files

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/param
file_temp="~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_ca"
declare -a Pd 
for i in {2..30..2}
do
        z=$i
        Pd+="$z " 
done
for m in {10..150..10}
do
for i in $Pd
do
        cd ~/ER_Ca_rot/networkFVMsims_git/examples/param_1D/len-$m/p-$i
        par=$(ls param*)

        for j in $par
        do
                timeout 10s ../../../../netmeshdynamicsFVM.exe "${j/param./}" 

        done
done
done