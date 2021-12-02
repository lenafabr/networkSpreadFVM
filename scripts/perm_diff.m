%perm difference plotter
addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");
dirname = "~/ER_Ca_rot/networkSpreadFVM-main/examples/inf_Perm_op/mesh_4";
dirname1 = "~/ER_Ca_rot/networkSpreadFVM-main/examples/200_perm_op";
ar = [0.0001,0.0005,0.0024,0.0111,0.0500];
ar = ar*0.3/0.0025;
s_a = string(ar);
l2 = length(ar);
dirn = strings(1,4);

dirinf = dir(dirname);
dirinf1 = dir(dirname1);
filen = {dirinf.name};
%filenn = {{'.'},{'..'},filen{4},filen{5},filen{6},filen{7},filen{8}};
filen1 = {dirinf1.name};

% dirinf = dir(dirname);
% filen = {dirinf.name};
% l = length(filen);

%%
scl = 1e-3*6e23/(1e3*1e12)/1e7;

cumflux = cell(2,l2);
time_avg = cell(2,l2);


    for j = 3:l2+2 %Kd
    j
    [fluxsim,tvals] = loadTotFluxSim(dirname+"/"+filen{j});
    if (tvals(1)~=0)
        tvals = [0; tvals];
        fluxsim = [fluxsim(1); fluxsim];
    end
    vals = sum(fluxsim,2);  
    integ = (vals(2:end)+vals(1:end-1))/2;
    tvals(1)
    dt = diff(tvals);
    cumflux0 = cumsum(integ.*dt);
    tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
   
    
    [fluxsim1,tvals1] = loadTotFluxSim(dirname1+"/"+filen1{j});
    if (tvals1(1)~=0)
        tvals1 = [0; tvals1];
        fluxsim1 = [fluxsim1(1); fluxsim1];
    end
    vals1 = sum(fluxsim1,2);  
    integ1 = (vals1(2:end)+vals1(1:end-1))/2;
    dt1 = diff(tvals1);
    cumflux1 = cumsum(integ1.*dt1);
    tavg1 = (tvals1(2:end)+tvals1(1:end-1))/2;
    
    
    
    cumflux{1}{j-2} = scl*cumflux0;
    time_avg{1}{j-2} = tavg0; 
    
    
    cumflux{2}{j-2} = scl*cumflux1;
    time_avg{2}{j-2} = tavg1; 
    
    
    end
%%


%%
hold off;
t = time_avg{1}{1}(1:10000);
%we get 25 edges in the radius of 3D0 of permnode(1D0 edges in hex network), so mM instantly loss is
%25 * 0.007854 = 0.19635mM = 0.19635 * 10^-3 * 6.023 * 10^23 /10^7 ~ 0.012
ad = 0.012;
for i = [1:l2]%Kd
    %if(i==3)
    subplot(2,3,i)
    hold off;
    sr = "Kd="+s_a(i)+"mM";
    ca = cumflux{1}{i}(1:10000) + ad;
    dif = cumflux{2}{i}(1:10000) - (cumflux{1}{i}(1:10000));
    dif1 = dif -ad;
    plot(t,dif,'--','DisplayName',"Diff");
    hold on;
    plot(t,dif1,'--','DisplayName',"Diff_shifted");
    plot(t,ca,'DisplayName',"Shifted");
    for j = 1:2 %inf and 200
        js = string(j); 
        plot(t,cumflux{j}{i}(1:10000),'DisplayName',js);
        hold on;
    end
    title(sr);
    xlabel("Time(s)");
    ylabel("Ca release(#ions/1e7)");
    legend;
    %end
end
