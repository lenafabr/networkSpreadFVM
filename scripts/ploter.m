addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");
dirname = "~/ER_Ca_rot/networkSpreadFVM-main/examples/Kd_op/Output";
%filen= strings(5,51);
ar = [0, 1.089, 2.178, 2.755, 4.356]; %diffusion constants
ar1 = [-13.815, -13.539, -13.262, -12.986, -12.710, -12.433, -12.157, -11.881, -11.604, -11.328, -11.052, -10.775, -10.499, -10.223, -9.947, -9.670, -9.394, -9.118, -8.841, -8.565, -8.289, -8.012, -7.736, -7.460, -7.183, -6.907, -6.631, -6.355, -6.078, -5.802, -5.526, -5.249, -4.973, -4.697, -4.420, -4.144, -3.868, -3.591, -3.315, -3.039, -2.763, -2.486, -2.210, -1.934, -1.657, -1.381, -1.105, -.828, -.552, -.276, 0]; %Kd exp
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
cnt = 1;
cumflux = zeros(54500,l-2);
time_avg = zeros(54500,l-2);

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
    cumflux(:,cnt) = scl*cumflux0;
    tavg0 = (tvals(2:end)+tvals(1:end-1))/2;
    %size(tavg0)
    time_avg(:,cnt) = tavg0;
    fluxsim0 = fluxsim; tvals0 = tvals;
    cnt=cnt+1
end


%%
%2725 steps per second
e = exp(1);
arr = 2.5 * e.^(ar1);

arr = arr/(pi*0.5*0.5);%gives Kd values in mM, a = 50nm
arr = arr(1:end)
length(arr)
conc_k = zeros(l2,5,l1); %index 1 --> diffusion consts; index 2 --> time; index 3 --> Kd 
%cnt =0;
for i = 1:l-2
    cnnt= 1;
    for j = [1,2,5,10,20]
        
            n = mod(i,l1)+1;
            m = floor(i/(l1+1));
            conc_k(m+1,cnnt,n) = cumflux(2725*j,i);
            
        cnnt = cnnt+1;
        end
end
%%
size(conc_k)
%%
c = zeros(1,l1);
arrp = gobjects(l2*5,1);
cnt =1;

hold off;
for i = 1:l2
    sp = subplot(2,3,i);
    cnnt =1;
    for j = [1,2,5,10,20]
     
        for k = 1:(l1)
            c(1,k) = conc_k(i,cnnt,k);
        end
        if(j == 1 || j==5 || j==20)
        st = "time="+j+"sec";
        arrp(cnt) = plot(arr,c, 'DisplayName',st);
        hold on;
        end
       
        cnt = cnt+1;
        cnnt = cnnt+1;
    end
        xlabel("Kd mM");
        ylabel("Cum. Ca release(mM)");
        title("D="+ar(i));
        legend(sp, [arrp(cnt-5);arrp(cnt-3);arrp(cnt-1)]);
end

%legend(arrp);





%%
% 
% a = -8.9872+0.2304*[0:26];
% a = e.^a;
% b = linspace(0.000125,0.05,27);
% e^-8.5624

%%

