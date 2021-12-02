%perm_mesh
addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");
dirname = "~/ER_Ca_rot/networkSpreadFVM-main/examples/bg_op";

ar = [0.0001,0.0005,0.0024,0.0111,0.0500];
ar = ar*0.3/0.0025;
s_a = string(ar);
l2 = length(ar);
dirn = strings(1,4);
filen = strings(1,l2);
dirinf = dir(dirname);
for i = 3:l2+2
filen(i-2) = dirinf(i).name;end

l = length(filen);
filen
%%
scl = 1e-3*6e23/(1e3*1e12)/1e7;

cumflux = cell(l2);
time_avg = cell(l2);

    for j = 3:l2+2 %Kd

    [fluxsim,tvals] = loadTotFluxSim(dirname+"/"+filen(i));
    if (tvals(1)~=0)
        tvals = [0; tvals];
        fluxsim = [fluxsim(1); fluxsim];
    end
    vals = sum(fluxsim,2);  
    integ = (vals(2:end)+vals(1:end-1))/2;
    dt = diff(tvals);
    cumflux0 = cumsum(integ.*dt);
    tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
   
    cumflux{j} = scl*cumflux0;
    time_avg{j} = tavg0; 
 
    end
%%


%%
hold off;
%t = time_avg{1}{1}(1:10000);

for i = [1:l2]%Kd
    %if(i==3)
    %subplot(2,3,i)
    sr = "Kd="+s_a(i)+"mM";
        js = string(j);
        plot(time_avg{i},cumflux{i},'DisplayName',js);
        hold on;
    title(sr);
    xlabel("Time(s)");
    ylabel("Ca release(mM)");
    legend;
   % end
end
