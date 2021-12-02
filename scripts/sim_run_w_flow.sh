#! /bin/bash
#this is run simulations on pc3 and store output files

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/
cp ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/log_param_wf/*  ~/ER_Ca_rot/networkSpreadFVM-main/examples/ 
par=$(ls ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/log_param_wf)

for i in $par
do
	timeout 10s ../netmeshdynamicsFVM.exe "${i/param./}" 

done
