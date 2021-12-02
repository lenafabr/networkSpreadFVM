addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
Kd=10.^[-2.045:0.2045:1.0225];
length(Kd)
t = 2.^[-5:1:10];
length(t)
p_c=20;
p=61;
y_e = zeros(16,16);
y0=3;
cnt1 = 1;
cnt2 = 1;
%%
for i = t
    cnt2=1;
    for j = Kd
       
        [x,y] = ode45(@(x,y) (-(p*y)/(1+((j*p_c)/(y+j)^2))) , [0,i], y0);
        y_e(cnt1,cnt2) = y(end);
        cnt2=cnt2+1;
        
    end
    cnt1=cnt1+1;
end
ind = zeros(length(t));
t_max = zeros(length(t));
c=1;
for i = [1:length(t)]
    length(y_e(i,:))
    [y_e2,ind(i)]=min(y_e(i,:));
    t_max(i) = Kd(ind(i));
    f=figure('visible','off');
    y_e1 = y_e(i,:);
    h = heatmap(y_e1);
    h.Title="time"+t(i);
    saveas(f,fullfile("~/ER_Ca_rot/networkSpreadFVM-main/examples/outputs","heat_map"+i), "jpeg");
    c=c+1
end
%%
length(t_max)
length(t)

f=figure('visible','on');
plot(log2(t),(t_max))
