#! /bin/bash
#this is run simulations on pc3 and store output files

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/param
file_temp="~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_ca"
declare -a Pd 
for i in {0..10}
do
	((z=(10+(2*i))))
	Pd+="$z " 
done

declare -a log_Kd
for i in {0..15}
do
	t=$(echo "scale=9;(2.045757491/10*$i)-2.045757491" | bc)
	log_Kd+="$t "
done

echo $log_Kd

cd ..

for i in $Pd
do 
	for j in $log_Kd
	do


		cp param.example1_ca param.to_change1
		po="STARTCONC 0.007854D0 0.15708D0"
		#echo $i
		p=$(echo "scale=9;(0.15708/20)*$i" | bc)
		#echo $j
		pn="STARTCONC 0.007854D0 ${p}D0"
		sed -i "s|$po|$pn|g" param.to_change1
		ko="KDEQUIL 2.5e-3"
		l_t=$(echo "scale=9;l(10)" | bc -l)
		k=$(echo "scale=9;8.33333*e($j*$l_t) " | bc -l)
		#lab=$(echo "scale=3;k/8.33" | bc)
		kn="KDEQUIL ${k}e-3"
		st="param.example1_ca_${i}_${k}e-3"
		ro="RESVINFO 1 5D3" 
       	rn="RESVINFO 1 5D3"
		sed -i "s|$ko|$kn|g" param.to_change1 
		sed -i "s|$ro|$rn|g" param.to_change1
		mv param.to_change1 ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/log_param_wof/bey-1/$st
done
done

