%perm_mesh
addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");
dirname = "~/ER_Ca_rot/networkSpreadFVM-main/examples/inf_Perm_op";

ar = [0.0001,0.0005,0.0024,0.0111,0.0500];
ar = ar*0.3/0.0025;
s_a = string(ar);
l2 = length(ar);
dirn = strings(1,4);
filen = cell(1,4);
for i = [1:4]
    num = string(i);
    dirn(i) = dirname+"/mesh_"+num;
    dirinf = dir(dirn(i));
    filen{i} = {dirinf.name};
end
% dirinf = dir(dirname);
% filen = {dirinf.name};
% l = length(filen);

%%
scl = 1e-3*6e23/(1e3*1e12)/1e7;

cumflux = cell(4,l2);
time_avg = cell(4,l2);

for i = 1:4 %mesh size
    for j = 3:l2+2 %Kd

    [fluxsim,tvals] = loadTotFluxSim(dirn(i)+"/"+filen{i}{j});
    if (tvals(1)~=0)
        tvals = [0; tvals];
        fluxsim = [fluxsim(1); fluxsim];
    end
    vals = sum(fluxsim,2);  
    integ = (vals(2:end)+vals(1:end-1))/2;
    dt = diff(tvals);
    cumflux0 = cumsum(integ.*dt);
    tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
   
    cumflux{i}{j-2} = scl*cumflux0;
    time_avg{i}{j-2} = tavg0; 
 
    end
end
%%


%%
hold off;
t = time_avg{1}{1}(1:10000);

for i = [1:l2]%Kd
    if(i==3)
    %subplot(2,3,i)
    hold off;
    sr = "Kd="+s_a(i)+"mM";
    for j = 1:4 %mesh
        js = string(j);
        plot(t,cumflux{j}{i}(1:10000),'DisplayName', "Mesh-0."+js);
        hold on;        
    end
    title(sr);
    xlabel("Time(s)");
    ylabel("Ca release(mM)");
    legend;
    end
end
