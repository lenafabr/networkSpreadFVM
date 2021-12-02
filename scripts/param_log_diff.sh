#! /bin/bash
#log_Kd and different protein diff
cd ~/ER_Ca_rot/networkFVMsims_git/examples/log_param_wof

declare -a prot
for i in {5..15}
do
	p=$(expr $i \* 2)
	prot+="$p "
done

declare -a diff
diff=("1" "3" "5")

for i in $prot
do
	cd p-$i
	parm_fl_edt=$(ls param.example1_ca_*)
	for j in $diff 
	do
	cd diff-$j
		for k in $parm_fl_edt
		do
			cp ../$k edit
        	ko="DCOEFF 27.55D0 2.755D0" 
        	kn="KDEQUIL 27.55D0 ${j}"

			st="${k}_diff_${j}"
        	sed -i "s|$ko|$kn|g" edit
        	mv edit $st
        	rm edit        	
		done
		cd .. #out of diff into p-
	done
	cd .. #out of p- into log_param
done




