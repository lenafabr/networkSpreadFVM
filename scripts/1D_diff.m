%1D diff all units in um,sec and #ions
num_cell=1000;
lth = 40;%um
lth_cel = lth/num_cell;%length per cell in um
rad=0.05;
vol_c=pi*rad^2 *(lth/num_cell); %um^3
cels=linspace(0,40,num_cell);
c_int_ca = 1;%mM
c_int_pb = 20;%mM
n_int_ca = c_int_ca*6.023*10^23 *vol_c*10^-18;
n_int_pb = c_int_pb*6.023*10^23 *vol_c*10^-18;
cel=zeros(num_cell,2);%each keeps track of unbound Ca and prot 
cel(:,1)=n_int_ca;
cel(:,2)=n_int_pb;
del_t = 0.001;%sec
% D_c = 2.8um^2/s, D_p=0.28um^2/s
t_u = 20;
t=0;
perm_cel = 0;
perm_rad = 3;%um
%%here we want the permeability to be infinite so the free ca in the
%%permeable cells is fixed to zero

while(t<=t_u)
    
    
    
    
    
    

    
t=t+del_t;
end