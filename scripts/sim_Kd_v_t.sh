#! /bin/bash
#this is run simulations on pc3 and store output files

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/
cp ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/Total_v_Kd_param/*  ~/ER_Ca_rot/networkSpreadFVM-main/examples/
par=$(ls ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/Total_v_Kd_param)

for i in $par
do
        timeout 10s ../netmeshdynamicsFVM.exe "${i/param./}" 

done
