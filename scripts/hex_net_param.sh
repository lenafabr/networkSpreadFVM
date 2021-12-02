#! /bin/bash

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/param
file_temp="~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_ca"
#declare -a Pd 
##do
#	((z=(10+(2*i))))
#	Pd+="$z " 
#done

#declare -a log_Kd
#for i in {0..10}
#do
#	t=$(echo "scale=9;((2.045757491/10)*$i)-2.045757491" | bc)
#	log_Kd+="$t "
#done

#echo $log_Kd
cd ..

for i in {1..6}
do 


		cp param.example1_ca param.to_change1
		po="NETFILE circlenuchexresv.net"
		#echo i
		pn="NETFILE hex_${i}_rd1.net"
		sed -i "s|$po|$pn|g" param.to_change1
		st="param.example1_hex_${i}_rd1"
		mv param.to_change1 ~/ER_Ca_rot/networkSpreadFVM-main/examples/$st
done
