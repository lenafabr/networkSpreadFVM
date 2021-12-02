#! /bin/bash
#this is run simulations on pc3 and store output files

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/
par=$(ls ~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_hex*)

for i in $par
do
	st="/home/aman/ER_Ca_rot/networkSpreadFVM-main/examples/param."
	timeout 10s ../netmeshdynamicsFVM.exe "${i/$st/}" 

done
