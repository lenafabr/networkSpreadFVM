%%
%dirarr = strings(54);
addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");
filearr = strings(54);
cnt = 1;
scl = 1e-3*6e23/(1e3*1e12)/1e7;
len = zeros(54,1);
cumflux = cell(11,11);
t_map = zeros(7,11,16);%time,prot,Kd
kd = [1:8];
kd = [1:8]*0.3/2.5;
K_dl = ((2.045757491/10)*[0:15])-2.045757491;
K_d = 8.3333*10.^(K_dl);
K_d = K_d * 0.3/2.5;
%cflx = cell(1,54);

%%
for i = 10:2:30
    n = (i/2)-4;
    d = "~/ER_Ca_rot/networkFVMsims_git/examples/heat_map_op/bey-1/p-"+i+"/";
    str = dir(d+"/*.out");
    dirn = {str.name};
    %get order right
    nnm=dirn(7:end);
    dirn{8}=nnm{3};
    dirn{9}=nnm{5};
    dirn{10}=nnm{7};
    dirn{11}=nnm{9};
    dirn{12}=nnm{2};
    dirn{13}=nnm{4};
    dirn{14}=nnm{6};
    dirn{15}=nnm{8};
    for j = 1:16
       fl = dirn{j};
        [fluxsim,tvals] = loadTotFluxSim(d+fl);
        if (tvals(1)~=0)
            tvals = [0; tvals]; 
            fluxsim = [fluxsim(1); fluxsim];
        end
        vals = sum(fluxsim,2);
        integ = (vals(2:end)+vals(1:end-1))/2;
        dt = diff(tvals);
        cumflux0 = cumsum(integ.*dt);
        cumflux{n,j} = cumflux0;
        %len(cnt) = (cumflux0(end));
        tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
        for k = 1:7
           dis = (abs(cumflux0-k)); 
           ind = (find(dis == min(dis)));
           t_map(k,n,j) = tavg0(ind);
           
            
        end
        
        fluxsim0 = fluxsim; tvals0 = tvals;
        %flux(((i/2)-4),j) = cumflux0(29509)*scl;
        
        fl = "~/ER_Ca_rot/networkSpreadFVM-main/examples/heat_map_op/";
        f = figure('visible','off');
        
        h0 = plot(tavg0,cumflux0*scl,'r-','LineWidth',3);
        xlabel('time (s)');
        ylabel('cumulative Ca ions released');
        s = "P-"+i+" Kd-"+K_d(j);
        title(s);
        xlim([0,20]);
        set(f,'visible','off');
        saveas(f,fullfile(fl,"time_plt_bey1"+s), "jpg");
        close(f);
        cnt=cnt+1
        
        
    end
        
        
        
  
end


cnt = 1;
%%

%%
for t = 1:7
yval = [10:2:30];
xval = K_d(5:end-2);
f1 = figure('visible','off');
h = heatmap(xval,yval(2:end),squeeze(t_map(t,2:end,5:end-2)));
axs = struct(gca);
cb = axs.Colorbar;
cb.Title.String = "sec";
h.Title = "Time for Ca o/p for "+t+"x10^7 ions";
h.YLabel = "Protein conc.(mM)";
h.XLabel = "Kd(mM)";
saveas(f1,fullfile(fll,"heat_map_bey1"+t),"jpg");
end
%%
flux;