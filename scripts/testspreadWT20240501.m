% load in and visualize snapshots
load('../../ERCaSims/networks/WT/WT_circreg.mat')
NT = WTnetworks(1);
%% load snapshot data

resdir = '../test/';
%filename = [resdir 'testrecoveryWT.mesh.txt'];
filename = [resdir 'testspreadWT.mesh.txt'];
MSH = MeshObj(filename);

%snapfile = [resdir 'testrecoveryWT.snap.txt'];
snapfile = [resdir 'testspreadWT.snap.txt'];
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
    %hold all    
%     cent =MSH.pos(centcell,:);
%     th=linspace(0,2*pi,100)';
%     circ = cent+Rabs*[cos(th) sin(th)];
%     plot(circ(:,1),circ(:,2),'g-','LineWidth',2)
    colormap copper
        
    caxis([0,0.1])

    %title(sprintf('Snap %d time %f', sc, snaptimes(sc)))
    % title('WT region 1')    
   
    set(gca,'Visible','off')    
    colorbar 
    plot_cleanup(gca,'FontSize',14,'pcolor',true)
    
    axis equal       
    drawnow
    
end
hold off

%% load concs and check total concentration over time
%snapfile = [resdir 'testrecoverynuc.snap.txt'];
snapfile = [resdir 'testspreadWT.snap.txt'];
[field,snaptimes,vels] = loadSnapshotFVM(snapfile);
nsnap = length(snaptimes)

% plot avg conc over time
scl = 1/(pi*0.05^2);
Ufield = squeeze(field(1:end,1,:));
Sfield = squeeze(field(1:end,2,:));
Cfield = Ufield.*(1 + Sfield./(Ufield+Kd));
%avgU = sum(Ufield.*MSH.len(1:end),1)/sum(MSH.len(1:end));

% get total calcium
Kd = 0.001583D0; % dissociation constant
S = 0.021342D0; % total binding sites
avgC = sum(Cfield.*MSH.len(1:end),1)/sum(MSH.len);
avgS = sum(Sfield.*MSH.len)/sum(MSH.len);
%avgC = avgU.*(1+S./(avgU + Kd));

% compare sims to analytical model, for both luminal and cytoplasmic calcium

a = 0.05;
Vtot = sum(MSH.len(1:end-1));
Vg = MSH.len(end);
Atot = sum(MSH.len(1:end-1));

ind = find(MSH.nodeind>0 | MSH.edgeind(:,1)'>0);
Anuc = 0;
Atot = sum(MSH.len(ind))+Anuc;
Vnuc = 0;
Vtot = sum(MSH.len(ind))+Vnuc;

plot(snaptimes, avgC*scl,'.-')
xlabel('time (s)')
ylabel('average ER total calcium')
legend off


%% make a linear network
nodepos = zeros(10,2);
nodepos(:,1) = 1:10;
edgenodes = [1:9 2:10];
NT = NetworkObj();
NT.nodepos = nodepos;
NT.edgenodes = edgenodes;
NT.setupNetwork();
NT.outputNetwork('../networks/linear10.net')