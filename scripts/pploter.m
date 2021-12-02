addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");
dirname = "~/ER_Ca_rot/networkSpreadFVM-main/examples/Kd_op/Op_new";
%filen= strings(5,51);
ar = [0, 0.275, 1.378, 2.755, 5.51, 27.55];%diffusion constants
st_ar = string(ar);
ar1 = [-8.9872:0.2304:-2.9968];%Kd exp
%ar1 = 2.5e-{element}
dirinf = dir(dirname);
filen = {dirinf.name};
l = length(filen);
l1 = length(ar1); %no of Kd values
l2 = length(ar); %no of diff constants

% for i = ar
%     count =1;
%     for j = ar1
%         %no = (-6+(j*3/25))/0.4343;
%         %round(no,4)
%         st = "example1_"+i+"_2.5e"+j+".out";
%         filen(find(ar==i),count) = dirname+st;
% end
% end


%%
scl = 1e-3*6e23/(1e3*1e12)/1e7;
ca = 27*[0:1:5];
cumflux = cell(l-2,1);
time_avg = cell(l-2,1);
cnt = 0;
for i = 3:length(filen)

    [fluxsim,tvals] = loadTotFluxSim(dirname+"/"+filen{i});
    st = dirname+"/"+filen{i};
    
    if (tvals(1)~=0)
        tvals = [0; tvals];
        fluxsim = [fluxsim(1); fluxsim];
    end
    vals = sum(fluxsim,2);  
    integ = (vals(2:end)+vals(1:end-1))/2;
    dt = diff(tvals);
    cumflux0 = cumsum(integ.*dt);
    %size(cumflux0)
    
    tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
    %size(tavg0)
    
    fluxsim0 = fluxsim; tvals0 = tvals;
    cnt=cnt+1
    
    
    if(contains(st, st_ar(2)))
        ca(2) = ca(2)+1;
        cumflux{ca(2)} = scl*cumflux0;
        time_avg{ca(2)} = tavg0;
    elseif(contains(st, st_ar(3)))
        ca(3) = ca(3)+1;
        cumflux{ca(3)} = scl*cumflux0;
        time_avg{ca(3)} = tavg0;
    elseif(contains(st, st_ar(4)))
        ca(4) = ca(4)+1;
        cumflux{ca(4)} = scl*cumflux0;
        time_avg{ca(4)} = tavg0;
    elseif(contains(st, st_ar(5)))
        ca(5) = ca(5)+1;
        cumflux{ca(5)} = scl*cumflux0;
        time_avg{ca(5)} = tavg0;
    elseif(contains(st, st_ar(6)))
        ca(6) = ca(6)+1;
        cumflux{ca(6)} = scl*cumflux0;
        time_avg{ca(6)} = tavg0;
    elseif(contains(st, st_ar(1)))
        ca(1) = ca(1)+1;
        cumflux{ca(1)} = scl*cumflux0;
        time_avg{ca(1)} = tavg0;
    end
            
    
end

%%

%%
%2725 steps per second
e = exp(1);
arr = e.^(ar1);

arr = arr/(pi*0.05*0.05);%gives Kd values in mM, a = 50nm


conc_k = zeros(l2,5,l1); %index 1 --> diffusion consts; index 2 --> time; index 3 --> Kd 
%cnt =0;

cnnt= 1;
for i = 1:l-2
    cnnt =1;
    %cflx = cumflux{i};
    %t = time_avg{i};
    for j = [1,2,5,10]
            n = mod(i,l1)+1;
            m = floor(i/(l1+1));
            cflx = cumflux{i};
            
            conc_k(m+1,cnnt,n) = cflx(2725*j);
     
         cnnt =cnnt+1;
    end

end

%%
size(conc_k);

% %%
% cumflux;
% c = zeros(1,l1);
% d = zeros(1,l1);
% for p = 1:l1
%    c(1,p) = conc_k(1,1,p); 
%    d(1,p) = conc_k(2,1,p); 
%     
% end
% plot(arr,c, arr,d);
%%
c = zeros(1,l1);
arrp = gobjects(l2*5,1);
cnt =1;
cp = 0;
hold off;
for i = 1:l2
    hold off;
    sp = subplot(2,3,i);
    cnnt =1;
    for j = [1,2,5,10]
     
        for k = 1:(l1)
            c(1,k) = conc_k(i,cnnt,k);
        end
        
        st = "time="+j+"sec";
        arrp(cnt) = plot(arr,c,'.','DisplayName',st);
        hold on;
        
       
        cnt = cnt+1;
        cnnt = cnnt+1;
    end
        xlabel("Kd mM");
        ylabel("Cum. Ca release(mM)");
        title("D="+ar(i));
        %legend(sp, [arrp(cnt-5);arrp(cnt-4);arrp(cnt-3);arrp(cnt-2);arrp(cnt-1)]);
        cp = cp +1
end
%legend(arrp);
% %%
% dirname = "~/ER_Ca_rot/networkSpreadFVM-main/examples/Kd_op/Op_new";
% 
% ff1  = dirname+"/example1_27.55_2.5e-5.0704.out";
% ff2 = dirname+"/example1_27.55_2.5e-5.9920.out";
% 
% filename = [ff1];
% f1 = [ff2];
% [fluxsim,tvals] = loadTotFluxSim(filename);
% [fl1,t1] = loadTotFluxSim(f1);
% if (tvals(1)~=0)
%     tvals = [0; tvals];
%     fluxsim = [fluxsim(1); fluxsim];
% end
% vals = sum(fluxsim,2);
% integ = (vals(2:end)+vals(1:end-1))/2;
% dt = diff(tvals);
% cumflux0 = cumsum(integ.*dt);
% tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
% fluxsim0 = fluxsim; tvals0 = tvals;
% 
% if (tvals(1)~=0)
%     t1 = [0; t1];
%     fl1 = [fl1(1); fl1];
% end
% vals = sum(fl1,2);
% integ = (vals(2:end)+vals(1:end-1))/2;
% dt = diff(t1);
% cumflux1 = cumsum(integ.*dt);
% tavg1 = (t1(2:end)+t1(1:end-1))/2;
% fluxsim1= fl1; tvals1 = t1;
% % plot cumulative flux
% scl = 1e-3*6e23/(1e3*1e12)/1e7;
% 
% h0 = plot(tavg0,cumflux0*scl,tavg1,cumflux1*scl,'LineWidth',3)
%%

