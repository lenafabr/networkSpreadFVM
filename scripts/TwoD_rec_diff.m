%2D circular 
%diffusion
n = 501;
r_l = linspace(0,50,n);% length scale--> 1/n 50um
r_s = 50/n; %um (0.1) seperation b/w cells
Uo = zeros(1,n);
Kd = 0.3;%mM
Du = 28; % diff constant L^2*T^-1 um^2/sec
Uo = Uo+1; %initial condition of unbound Ca as 1mM
Un = Uo;
Po = zeros(1,n);
Po = Po+50;%inittial condition of protein buffer as 50mM
Pn = Po;
Dp = 2.8;%diff constant of buffer um^2/s
%so relevant time scale --> (D*(1/n)^-2)^-1; but timestep can be arbitary
t_s = ((Du*((r_s)^-2))^-1)/10; %div by 10 to make sure of stability ~(1/2800 sec * 1/10)
t_stps = floor(1/t_s);
st = "time_step" + 0;
%plot(r_l,Uo,r_l,Po,'DisplayName',st);
%xlabel("r um");
%ylabel("Conc mM");
J = zeros(1,t_stps);
no_ca = zeros(1,t_stps);
cflx = zeros(1,t_stps);
%%
t_stps;
%%

for i = 1:t_stps %time steps 28k for 1 sec 
    for j = 1:n %going over space 0.1 um cell at a time
        if(j==1)
            Un(j) = 0; %say absorbing boundaries
            Pn(j) = 0;
        elseif(j~=n)
            delUn =((1/(1+(Pn(j)*Kd/(Un(j)+Kd)^2)))*((Un(j+1)-2*Un(j)+Un(j-1))/r_s^2 * (Du + (Dp*Kd*Pn(j))/(Un(j)+Kd)^2) + ((Un(j+1)-Un(j))/r_s)*(1/(j*r_s))*((Dp*Kd*Pn(j))/(Un(j)+Kd)^2) - 2*((Un(j+1)-Un(j))/r_s)*((Dp*Kd*Pn(j))/(Un(j)+Kd)^3)))*t_s;
            delPn = (Dp*(Pn(j+1)-2*Pn(j)+Pn(j-1))/r_s^2 + Dp/(r_s*j) * (Pn(j+1)-Pn(j))/r_s)*t_s;
            Un(j) = Un(j) +delUn;
            Pn(j) = Pn(j) +delPn; 
            %if Cn(j)<0 ; Cn(j) =0;end
        end
    end
    J(i) = Du*(Un(2)-Un(1))/r_s;
    no_ca(i) = J(i)*2*pi*r_s;
    cflx(i) = sum(no_ca);
%     figure('visible','off');
%     if(mod(i,1000) ==0)
%         t = r_s * j;
%         stC = "time_" +t+"sec_Ca";
%         p1 = plot(r_l,Un,'DisplayName',stC);
%         saveas(p1,stC+".png");
%         stP = "time_" +t+"sec_Pb";
%         p2 = plot(r_l,Pn,'DisplayName',stP);nano
%         saveas(p2,stP+".png");
%     end
    i
end
p3 = plot(1:t_stps,cmflx);
