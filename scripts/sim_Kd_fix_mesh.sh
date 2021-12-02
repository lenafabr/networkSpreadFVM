#! /bin/bash
#this is run simulations on pc3 and store output files
for j in {1..4}
do
cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/
cp ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/perm_fix/mesh_$j/*  ~/ER_Ca_rot/networkSpreadFVM-main/examples/ 
par=$(ls ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/perm_fix/mesh_$j)

for i in $par
do
	timeout 335s ../netmeshdynamicsFVM.exe "${i/param./}" 

done
done