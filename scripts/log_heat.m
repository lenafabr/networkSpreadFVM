%%
%log_heat
%To get cumalative Ca flux by varying the proteins conc. and Kd

addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
filearr_lwf = strings(11,11);
filearr_lwof=strings(11,11);
fa_lwof = strings(11);
fa_lwf = strings(11);
cnt = 1;
scl = 1e-3*6e23/(1e3*1e12)/1e7;
len_wf = zeros(121);
len_wof=zeros(121);
flux_wf= zeros(11,11);
flux_wof=zeros(11,11);

%%
for i = 10:2:30
    str = "~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wf_op/P-"+i+"mM/";
    str1 = "~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wof_op1/log_wof_op/P-"+i+"mM/";
    str2 = "~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wof_op1/log_wof_op/P-"+i+"mM/";
    %cd(str);
    
   % inf=(dir('*.out'));
   % filearr_lwf(cnt,:)={inf.name};
    cd(str1);
    inf1=(dir('*.out'));
    filearr_lwof(cnt,:)={inf1.name};
 
   % cd(str2);
   % inf2=(dir('*.out'));
   % fa_lwof={inf2.name};
        
        

       %{ 
        
        f = figure('visible','off');
        
        h0 = plot(tavg0,cumflux0*scl,'r-','LineWidth',3);
        xlabel('time (s)');
        ylabel('cumulative Ca ions released');
        s = "P-"+i+"Kd-"+j;
        title(s);
        xlim([0,20]);
        saveas(f,fullfile(fl,"time_plt"), "fig");
        close(f);
         %}
    
   cnt=cnt+1;
end

%%
filearr_lwof
%%
cnt=1;
for i = 10:2:30
    %str = "~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wf_op/P-"+i+"mM/";
    str1 = "~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wof_op1/log_wof_op/P-"+i+"mM/";
    for j = 1:11
        
        n = ((i-10)/2)+1;
        %{
        [fluxsim,tvals] = loadTotFluxSim(str+filearr_lwf(n,j));
        if (tvals(1)~=0)
            tvals = [0; tvals]; 
            fluxsim = [fluxsim(1); fluxsim];
        end
        vals = sum(fluxsim,2);
        integ = (vals(2:end)+vals(1:end-1))/2;
        dt = diff(tvals);
        cumflux0 = cumsum(integ.*dt);
       
        len_wf(cnt)=length(cumflux0);
        flux_wf(n,j)=cumflux0(22010)*scl;
        %}
        

        [fluxsim,tvals] = loadTotFluxSim(str1+filearr_lwof(n,j));
        if (tvals(1)~=0)
            tvals = [0; tvals]; 
            fluxsim = [fluxsim(1); fluxsim];
        end
        vals = sum(fluxsim,2);
        integ = (vals(2:end)+vals(1:end-1))/2;
        dt = diff(tvals);
        cumflux0 = cumsum(integ.*dt);
        len_wof(cnt)=length(cumflux0);
        flux_wof(n,j)=cumflux0(27761)*scl;
        cnt=cnt+1
    end
end

%%
flux=zeros(11);
l=zeros(11);
length(fa_lwof)
cnt=1;
for j = 1:11
    [fluxsim,tvals] = loadTotFluxSim(str2+fa_lwof(j));
        if (tvals(1)~=0)
            tvals = [0; tvals]; 
            fluxsim = [fluxsim(1); fluxsim];
        end
        vals = sum(fluxsim,2);
        integ = (vals(2:end)+vals(1:end-1))/2;
        dt = diff(tvals);
        cumflux0 = cumsum(integ.*dt);
        l(j)=length(cumflux0);
        flux(j)=cumflux0(25802)*scl;
end
min(l)
%%
x = [-2.045:0.2045:0];
xval = 10.^x;
yval = [10:2:30];
f = figure('visible','on');
h = heatmap(xval,yval,flux_wof);
h.Title = "Without flow";
h.XLabel ="K_d(mM)";
h.YLabel = "Prot conc.(mM)";
%%


cnt=1;
xval = [-1.045:0.1045:0];
yval = [10:2:30];

        
f = figure('visible','off');
h = heatmap(xval,yval,flux_wf);
h.Title = "With Flow";
h.YLabel = "Protein conc.(mM)";
h.XLabel = "log_Kd";
saveas(f,fullfile("~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wf_op/","heat_map"), "jpeg");
close(f);
        
f = figure('visible','off');
h = heatmap(xval,yval,flux_wof);
h.Title = "Without flow";
h.YLabel = "Protein conc.(mM)";
h.XLabel = "log_Kd";
saveas(f,fullfile("~/ER_Ca_rot/networkSpreadFVM-main/examples/log_wof_op/","heat_map"), "jpeg");
close(f);     

%min(len_wf)
%min(len_wof)
%26561 wof
%22010 wf

        %{
        tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
        fluxsim0 = fluxsim; tvals0 = tvals;
        flux(((i/2)-4),j) = cumflux0(29509)*scl;%}


        
%plotting heat map
xval = [10:2:22];
yval = [1:9];
f = figure('visible','on');
h = heatmap(yval,xval,flux);
h.Title = "Final cumilative flux output for Ca";
h.YLabel = "Protein conc.(mM)";
h.XLabel = "Kd(mM)";
%%
csvwrite('~/ER_Ca_rot/networkSpreadFVM-main/examples/heat_map_op/heat_map.txt',flux);
%%
%plot Ca vs Kd/[P]
ct =1;
for i = 1:7
    for j = 1:9
        
        x = 1:7;
        x1 = 1:9;
        
        f = figure('visible','off');
        scatter(x1,flux(i,:));
        fl = "~/ER_Ca_rot/networkSpreadFVM-main/examples/heat_map_op/";
        s = "Kd"+i;
        saveas(f,fullfile(fl,s), "jpeg");
        close(f);
        f = figure('visible','off');
        scatter(x,flux(:,j));
        fl = "~/ER_Ca_rot/networkSpreadFVM-main/examples/heat_map_op/";
        s = "[P]"+j;
        saveas(f,fullfile(fl,s), "jpeg");
        close(f);
        ct = ct+1
    end
end

        %}