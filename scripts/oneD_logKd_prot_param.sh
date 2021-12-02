cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/
file_temp="~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_ca"
declare -a Pd 
for i in {2..30..2}
do
	z=$i
	Pd+="$z " 
done

declare -a log_Kd
for i in {0..12}
do
	t=$(echo "scale=9;(1.8/10*$i)-1.8" | bc)
	log_Kd+="$t "
done

echo $log_Kd

cd ..
for m in {10..150..10}
do
for i in $Pd
do 
	for j in $log_Kd
	do


		cp param.exampl1_1D param.to_change1
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
		lo="NETFILE line.net"
		ln="NETFILE line${m}.net"
		sed -i "s|$ko|$kn|g" param.to_change1
		sed -i "s|$lo|$ln|g" param.to_change1
		mv param.to_change1 ~/ER_Ca_rot/networkFVMsims_git/examples/param_1D/len-$m/p-$i/$st
done
done
done
