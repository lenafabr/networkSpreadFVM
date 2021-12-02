#! /bin/bash

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/param
file_temp="~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_ca"
declare -a Pr 
for i in {0..20}
do
	((z=($i*10)))
	Pr+="$z " 
done


cd ..

for i in $Pr
do 


		cp param.example1_ca param.to_change1
		po="PERMNODE 1 200D0"
		echo $i
		#p=$(echo "scale=9;(0.15708/20)*$i" | bc)
		pn="PERMNODE 1 $i"
		sed -i "s|$po|$pn|g" param.to_change1
		st="param.example1_ca_perm_$i"
		mv param.to_change1 ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/perm/$st
done

