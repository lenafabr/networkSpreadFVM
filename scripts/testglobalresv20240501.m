% load in and visualize snapshots
%load('../../ERCaSims/networks/WT/WT_circreg.mat')
%NT = WTnetworks(1);
%
NT = NetworkObj('../test/circlenuchexresv.net',struct('dim',2));
NT.setCumEdgeLen(1:NT.nedge);
%% load snapshot data

resdir = '../test/';
%filename = [resdir 'testrecoveryWT.mesh.txt'];
filename = [resdir 'testrecoverynuc.mesh.txt'];
%filename = [resdir 'test.mesh.txt'];
MSH = MeshObj(filename);

%snapfile = [resdir 'testrecoveryWT.snap.txt'];
snapfile = [resdir 'testrecoverynuc.snap.txt'];
%snapfile = [resdir 'test.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes)

%[~,centcell] = min(field(:,1,2))
% TODO:
% FIGURE out why different mesh cells are getting different concentrations
% with global recovery.
% figure out why total ER concentration is not asymptoting correctly

%% show the snapshot
showtime = 2; % time at which to show in sec
 [~,showsnap] = min(abs(snaptimes - showtime));
 Rabs = 0.25;
 tubeR = 0.05;
 clf
for sc = 1:nsnap
    [sc nsnap]
    ind = find(MSH.resvind==0);
    ind(ind==MSH.globalresvind) = [];
    %intpos = interpolateMeshPos(MSH,NT);    
    %scatter(intpos(ind,1),intpos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
    scatter(MSH.pos(ind,1),MSH.pos(ind,2),15,field(ind,1,sc)/0.007854,'filled')   
    
    %[X,Y,patchind] = meshPatches(MSH,NT,tubeR,tubeR);
    %C = field(patchind,1,sc)/(pi*tubeR^2); % colors
    %patch(X,Y,C);
    %shading flat    
    hold all    
%     cent =MSH.pos(centcell,:);
%     th=linspace(0,2*pi,100)';
%     circ = cent+Rabs*[cos(th) sin(th)];
%     plot(circ(:,1),circ(:,2),'g-','LineWidth',2)
    colormap copper
        
    %caxis([0,0.5])
    cmax = max(field(ind,1,sc)/0.007854);
    if (abs(cmax)<eps); cmax = 0.1; end

%     % show nucleus
    nucind = find(MSH.resvind>0); 
    nuccent = MSH.pos(nucind,:);
    th = linspace(0,2*pi)';
    nucpts = nuccent + 7.5*[cos(th) sin(th)];           
    cmap = colormap(gca);
    nuccol = interp1(linspace(0,cmax,size(cmap,1)),cmap,field(nucind,1,sc)/(pi*0.05^2),'linear','extrap');
    fill(nucpts(:,1),nucpts(:,2),nuccol,'LineStyle','none','FaceAlpha',1)

    %title(sprintf('Snap %d time %f', sc, snaptimes(sc)))
    % title('WT region 1')    
   
    set(gca,'Visible','off')    
    colorbar 
    plot_cleanup(gca,'FontSize',14,'pcolor',true)
    
    axis equal       
    drawnow
    
end
hold off



%% load concs and compare to analytic results
snapfile = [resdir 'testrecoverynuc.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes)

% plot average conc over time
scl = 1/(pi*0.05^2);
Ufield = squeeze(field(1:end-1,1,:));
Sfield = squeeze(field(1:end-1,2,:));
Cfield = Ufield.*(1 + Sfield./(Ufield+Kd));
avgC = sum(Cfield.*MSH.len(1:end-1),1)/sum(MSH.len(1:end-1));
avgS = sum(Sfield.*MSH.len(1:end-1))/sum(MSH.len(1:end-1));
%avgU = sum(Ufield.*MSH.len(1:end-1),1)/sum(MSH.len(1:end-1));

% get total calcium
Kd =  0.001583D0; % dissociation constant
S = 0.021342D0; % total binding sites
%S = 0
%avgC = avgU.*(1+S./(avgU + Kd));

% compare sims to analytical model, for both luminal and cytoplasmic calcium

a = 0.05;
Vtot = sum(MSH.len(1:end-1));
Vg = MSH.len(end);
Atot = sum(MSH.len(1:end-1));

ind = find(MSH.nodeind>0 | MSH.edgeind(:,1)'>0);
Anuc = 2.25D3;
%Anuc = 0;
Atot = sum(MSH.len(ind))+Anuc;
Vnuc = 5.4D3;
%Vnuc = 0;
Vtot = sum(MSH.len(ind))+Vnuc;

cg0 = 6.4093D-4;

cfin = cg0*Vg/Vtot*scl;

kr = 100;
kout = 0;
ratrecover = 1;

subplot(1,2,1)
plot(snaptimes, avgC*scl,'.',snaptimes,ratrecover*field(end,1,1)*Vg/Vtot.*(1-exp(-(kr*Atot+kout)/Vg*snaptimes))*scl,'r-','LineWidth',2,'MarkerSize',15)
plot_cleanup(gca,'FontSize',14)
xlabel('time (s)')
ylabel('average ER total calcium')
legend off
subplot(1,2,2)
plot(snaptimes, squeeze(field(end,1,:))*scl,'b.',snaptimes,field(end,1,1).*exp(-(kr*Atot+kout)/Vg*snaptimes)*scl,'r-','LineWidth',2,'MarkerSize',15)
plot_cleanup(gca,'FontSize',14)
xlabel('time (s)')
ylabel('cytoplasmic calcium')
legend off