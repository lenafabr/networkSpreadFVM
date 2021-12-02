#! /bin/bash

cd ~/ER_Ca_rot/networkSpreadFVM-main/examples/param
file_temp="~/ER_Ca_rot/networkSpreadFVM-main/examples/param.example1_ca_fix"

declare -a mesh
for i in {1..4}
do
	p=$(echo "scale=4;($i*0.1)" | bc)
	mesh+="$p "
done



declare -a log_Kd 
for i in {0..4}
do
        z=$(echo "scale=4;-8.9872+($i*1.498)" | bc)

        log_Kd+="$z " 
done
#ts=("4" "10" "19" "38")
cd ..



 	
for i in $log_Kd
do
for j in $mesh
do
        cp param.example1_ca_fix param.to_change1
        ko="KDEQUIL 2.5e-3"
        #l_t=$(echo "scale=9;l(10)" | bc -l)
        #k=$(echo "scale=9;8.33333*e($j*$l_t) " | bc -l)
        #lab=$(echo "scale=3;k/8.33" | bc)
        num=$(echo "scale=4;e($i)" | bc -l) 
        kn="KDEQUIL $num"
	mo="MESHSIZE 0.2D0"
	mn="MESHSIZE ${j}D0"
	jn=$(echo "scale=0;$j*10" | bc -l)
	jn=${jn%.*}
	st="mesh_${jn}/param.example1_ca_m${jn}_perm_$num"
        sed -i "s|$ko|$kn|g" param.to_change1
	sed -i "s|$mo|$mn|g" param.to_change1
        mv param.to_change1 ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/perm_fix/$st


done
done


#cp param.example1_ca param.to_change
#ko="KDEQUIL 2.5e-3"
#kn="KDEQUIL 0"
#st="param.example1_2.5e-{i}"
#sed -i "s|$ko|$kn|g" param.to_change1 
#mv param.to_change1 ~/ER_Ca_rot/networkSpreadFVM-main/examples/param/Total_v_Kd_param
