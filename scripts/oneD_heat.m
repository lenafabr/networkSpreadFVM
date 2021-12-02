%%
%dirarr = strings(54);
addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");

cnt = 1;
scl = 1e-3*6e23/(1e3*1e12)/1e7;%0.06
len = zeros(54,1);
cumflux = cell(11,13);
t_map = zeros(7,11,13);%time,prot,Kd
kd = [1:8];
kd = [1:8]*0.3/2.5;
K_dl = ((1.8/10)*[0:12])-1.8;
K_d = 8.3333*10.^(K_dl);
K_d = K_d * 0.3/2.5;
%cflx = cell(1,54);

%%
for i = 10:2:30
    n = (i/2)-4;
    d = "~/ER_Ca_rot/networkFVMsims_git/examples/oneD_op/p-"+i+"/";
    str = dir(d+"/*.out");
    dirn = {str.name};

    %get order right
    nnm=dirn(8:end);
    dirn{8}=nnm{3};
    dirn{9}=nnm{4};
    dirn{10}=nnm{5};
    dirn{11}=nnm{6};
    dirn{12}=nnm{1};
    dirn{13}=nnm{2};

    for j = 1:length(dirn)
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
        f = cumflux0*scl;
        cumflux{n,j} = f;
        %len(cnt) = (cumflux0(end));
        tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
        for k = 1:7
           kn = 2^(-k/2); 
           dis = (abs(cumflux0-kn)); 
           ind = (find(dis == min(dis)));
           t_map(k,n,j) = tavg0(ind);
           
            
        end
      
        fluxsim0 = fluxsim; tvals0 = tvals;
        %flux(((i/2)-4),j) = cumflux0(29509)*scl;
        
        fl = "~/ER_Ca_rot/networkFVMsims_git/examples/oneD_op/plots/time_plts/";
        f = figure('visible','off');
        
        h0 = plot(tavg0,cumflux0*scl,'r-','LineWidth',3);
        xlabel('time (s)');
        ylabel('cumulative Ca ions released');
        s = "P-"+i+" Kd-"+K_d(j);
        title(s);
        xlim([0,20]);
        set(f,'visible','off');
        saveas(f,fullfile(fl,"time_plt"+s), "jpg");
        close(f);
        cnt=cnt+1
        
        
    end
        
        
        
  
end


cnt = 1;
%%
fll= "~/ER_Ca_rot/networkFVMsims_git/examples/oneD_op/plots/";
ll=length(tavg0);
for i = 10:2:30
     nm = (i/2)-4;
     arr=zeros(1,length(dirn));
    for j = 1:length(dirn)
        arr(j)=cumflux{nm,j}(ll);        
    end
    f=figure('visible','off');
    plot(K_d,arr);
    title("P-"+i+"mM");
    hold on;
    saveas(f,fullfile(fll,"KD-p-"+i+"mM"),"jpg");
end
%%
for t = 1:7
tn = 2^(-t/2)*scl;
yval = [10:2:30];
xval = K_d(1:end);
f1 = figure('visible','on');
h = heatmap(xval,yval(1:end),squeeze(t_map(t,1:end,1:end)));
axs = struct(gca);
cb = axs.Colorbar;
cb.Title.String = "sec";
h.Title = "Time for Ca o/p for "++"x10^7 ions";
h.YLabel = "Protein conc.(mM)";
h.XLabel = "Kd(mM)";
fll= "~/ER_Ca_rot/networkFVMsims_git/examples/oneD_op/plots/";
saveas(f1,fullfile(fll,"heat_map"+t),"jpg");
end