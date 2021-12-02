#! /bin/bash
#log_Kd and different protein diff
cd ~/ER_Ca_rot/networkFVMsims_git/examples/heat_map_op

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
	for j in $diff 
	do
	cd diff-$j
	parm_fl=$(ls ~/ER_Ca_rot/networkFVMsims_git/examples/log_param_wof/p-${i}/diff-${j}/param.example1_ca_*)
		for k in $parm_fl
		do
		s="~/ER_Ca_rot/networkFVMsims_git/examples/log_param_wof/p-${i}/diff-${j}/$k"	
		timeout 300s ../../../../../netmeshdynamics.exe "${s/param./}" 	
             	
		done
		cd .. #out of diff into p-
	done
	cd .. #out of p- into log_param
done




